/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2025 NKI/AVL, Netherlands Cancer Institute
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "revision.hpp"
#include "tortoize.hpp"

#include <zeep/crypto.hpp>
#include <zeep/http/daemon.hpp>
#include <zeep/http/html-controller.hpp>
#include <zeep/http/rest-controller.hpp>
#include <zeep/http/server.hpp>

#include <cif++/cif++.hpp>
#include <mcfp/mcfp.hpp>

#include <fstream>

namespace fs = std::filesystem;

using json = nlohmann::json;

#ifdef _MSC_VER
# include <fcntl.h>
// MSVC stdlib.h definitions
# define STDIN_FILENO 0
# define STDOUT_FILENO 1
# define STDERR_FILENO 2
#endif

// --------------------------------------------------------------------

class tortoize_html_controller : public zeep::http::html_controller
{
  public:
	tortoize_html_controller()
		: zeep::http::html_controller("tortoize")
	{
		mount("{css,scripts,fonts,images,favicon}/", &tortoize_html_controller::handle_file);
		mount("{favicon.ico,browserconfig.xml,manifest.json}", &tortoize_html_controller::handle_file);
		mount("", &tortoize_html_controller::index);
	}

	void index(const zeep::http::request &request, const zeep::http::scope &scope, zeep::http::reply &reply)
	{
		get_template_processor().create_reply_from_template("index", scope, reply);
	}
};

// --------------------------------------------------------------------

class tortoize_rest_controller : public zeep::http::rest_controller
{
  public:
	tortoize_rest_controller()
		: zeep::http::rest_controller("")
		, m_tempdir(fs::temp_directory_path() / "tortoize-ws")
	{
		fs::create_directories(m_tempdir);

		map_post_request("tortoize", &tortoize_rest_controller::calculate, "data", "dict");
	}

	json calculate(const std::string &file, const std::string &dict)
	{
		// First store dictionary, just in case

		fs::path dictFile;

		if (not dict.empty())
		{
			dictFile = m_tempdir / ("dict-" + std::to_string(m_next_dict_nr++));
			std::ofstream tmpFile(dictFile);
			tmpFile << dict;

			cif::compound_factory::instance().push_dictionary(dictFile);
		}

		try
		{
			// --------------------------------------------------------------------

			json data{
				{ "software",
					{ { "name", "tortoize" },
						{ "version", kVersionNumber },
						{ "reference", "Sobolev et al. A Global Ramachandran Score Identifies Protein Structures with Unlikely Stereochemistry, Structure (2020)" },
						{ "reference-doi", "https://doi.org/10.1016/j.str.2020.08.005" } } }
			};

			// --------------------------------------------------------------------

			struct membuf : public std::streambuf
			{
				membuf(char *text, size_t length)
				{
					this->setg(text, text, text + length);
				}
			} buffer(const_cast<char *>(file.data()), file.length());

			cif::gzio::istream in(&buffer);

			cif::file f = cif::pdb::read(in);
			if (f.empty())
				throw std::runtime_error("Invalid mmCIF or PDB file");

			std::set<uint32_t> models;
			for (auto r : f.front()["atom_site"])
			{
				if (not r["pdbx_PDB_model_num"].empty())
					models.insert(r["pdbx_PDB_model_num"].as<uint32_t>());
			}

			if (models.empty())
				models.insert(0);

			for (auto model : models)
			{
				cif::mm::structure structure(f, model);
				data["model"][std::to_string(model)] = calculateZScores(structure);
			}

			if (not dictFile.empty())
			{
				cif::compound_factory::instance().pop_dictionary();
				fs::remove(dictFile);
			}

			return data;
		}
		catch (...)
		{
			std::error_code ec;

			if (not dictFile.empty())
			{
				cif::compound_factory::instance().pop_dictionary();
				fs::remove(dictFile);
			}

			throw;
		}
	}

	fs::path m_tempdir;
	size_t m_next_dict_nr = 1;
};

int start_server(int argc, char *argv[])
{
	using namespace std::literals;
	namespace zh = zeep::http;

	cif::compound_factory::init(true);

	int result = 0;

	auto &config = mcfp::config::instance();

	config.init("tortoize server [options] start|stop|status|reload",
		mcfp::make_option("help,h", "Display help message"),
		mcfp::make_option("version", "Print version"),
		mcfp::make_option("verbose,v", "verbose output"),

		mcfp::make_option<std::string>("address", "0.0.0.0", "External address"),
		mcfp::make_option<uint16_t>("port", 10350, "Port to listen to"),
		mcfp::make_option<std::string>("user,u", "www-data", "User to run the daemon"),

		mcfp::make_option("no-daemon,F", "Do not fork into background"));

	config.parse(argc, argv);

	// --------------------------------------------------------------------

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.has("help"))
	{
		std::cout << config << std::endl;
		exit(0);
	}

	if (config.operands().empty())
	{
		std::cerr << "Missing command, should be one of start, stop, status or reload" << std::endl;
		exit(1);
	}

	cif::VERBOSE = config.count("verbose");

	std::string user = config.get<std::string>("user");
	std::string address = config.get<std::string>("address");
	uint16_t port = config.get<uint16_t>("port");

	zh::daemon server([&]()
		{
		 auto s = new zeep::http::server();

#if DEBUG
		 s->set_template_processor(new zeep::http::file_based_html_template_processor("docroot"));
#else
		 s->set_template_processor(new zeep::http::rsrc_based_html_template_processor());
#endif
		 s->add_controller(new tortoize_rest_controller());
		 s->add_controller(new tortoize_html_controller());
		 return s; }, "tortoized");

	std::string command = config.operands().front();

	if (command == "start")
	{
		std::cout << "starting server at http://" << address << ':' << port << '/' << std::endl;

		if (config.has("no-daemon"))
			result = server.run_foreground(address, port);
		else
			result = server.start(address, port, 2, 2, user);
		// server.start(vm.count("no-daemon"), address, port, 2, user);
		// // result = daemon::start(vm.count("no-daemon"), port, user);
	}
	else if (command == "stop")
		result = server.stop();
	else if (command == "status")
		result = server.status();
	else if (command == "reload")
		result = server.reload();
	else
	{
		std::cerr << "Invalid command" << std::endl;
		result = 1;
	}

	return result;
}

// --------------------------------------------------------------------

// recursively print exception whats:
void print_what(const std::exception &e)
{
	std::cerr << e.what() << std::endl;
	try
	{
		std::rethrow_if_nested(e);
	}
	catch (const std::exception &nested)
	{
		std::cerr << " >> ";
		print_what(nested);
	}
}

int main(int argc, char *argv[])
{
	int result = -1;

	try
	{
		result = start_server(argc, argv);
	}
	catch (std::exception &ex)
	{
		print_what(ex);
		exit(1);
	}

	return result;
}
