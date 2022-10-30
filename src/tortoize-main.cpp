/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
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

#include <fstream>

#include <zeep/http/daemon.hpp>
#include <zeep/http/server.hpp>
#include <zeep/http/html-controller.hpp>
#include <zeep/http/rest-controller.hpp>
#include <zeep/crypto.hpp>

#include <gxrio.hpp>
#include <cfg.hpp>
#include <cif++.hpp>

#include "tortoize.hpp"
#include "revision.hpp"

namespace fs = std::filesystem;
namespace ba = boost::algorithm;

using json = zeep::json::element;


#ifdef _MSC_VER
#include <fcntl.h>
//MSVC stdlib.h definitions
#define STDIN_FILENO 0
#define STDOUT_FILENO 1
#define STDERR_FILENO 2
#endif



// --------------------------------------------------------------------
#if WEBSERVICE

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

	void index(const zeep::http::request& request, const zeep::http::scope& scope, zeep::http::reply& reply)
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

	json calculate(const std::string& file, const std::string& dict)
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
					{
						{ "name", "tortoize" },
						{ "version", kVersionNumber },
						{ "reference", "Sobolev et al. A Global Ramachandran Score Identifies Protein Structures with Unlikely Stereochemistry, Structure (2020)" },
						{ "reference-doi", "https://doi.org/10.1016/j.str.2020.08.005" }
					}
				}
			};

			// --------------------------------------------------------------------

			struct membuf : public std::streambuf
			{
				membuf(char *text, size_t length)
				{
					this->setg(text, text, text + length);
				}
			} buffer(const_cast<char *>(file.data()), file.length());

			gxrio::istream in(&buffer);

			cif::file f = cif::pdb::read(in);
			if (f.empty())
				throw std::runtime_error("Invalid mmCIF or PDB file");

			std::set<uint32_t> models;
			for (auto r: f.front()["atom_site"])
			{
				if (not r["pdbx_PDB_model_num"].empty())
					models.insert(r["pdbx_PDB_model_num"].as<uint32_t>());
			}

			if (models.empty())
				models.insert(0);

			for (auto model: models)
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

int start_server(int argc, char* argv[])
{
	using namespace std::literals;
	namespace zh = zeep::http;

	cif::compound_factory::init(true);

	int result = 0;

	auto &config = cfg::config::instance();

	config.init("tortoize server [options] start|stop|status|reload",
		cfg::make_option("help,h", "Display help message"),
		cfg::make_option("version", "Print version"),
		cfg::make_option("verbose,v", "verbose output"),

		cfg::make_option<std::string>("address", "0.0.0.0", "External address"),
		cfg::make_option<uint16_t>("port", 10350, "Port to listen to"),
		cfg::make_option<std::string>("user,u", "www-data", "User to run the daemon"),

		cfg::make_option("no-daemon,F", "Do not fork into background"));

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
		return s;
	}, kProjectName );

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
#endif

// --------------------------------------------------------------------

int pr_main(int argc, char* argv[])
{
	using namespace std::literals;

#if WEBSERVICE
	if (argc > 2 and argv[1] == "server"s)
		return start_server(argc - 1, argv + 1);
#endif

	auto &config = cfg::config::instance();

	config.init("tortoize [options] input [output]",
		cfg::make_option("help,h", "Display help message"),
		cfg::make_option("version", "Print version"),

		cfg::make_option("verbose,v", "verbose output"),

		cfg::make_option<std::string>("log", "Write log to this file"),

		cfg::make_option<std::vector<std::string>>("dict",
			"Dictionary file containing restraints for residues in this specific target, can be specified multiple times."),

		cfg::make_hidden_option<std::string>("build", "Build a binary data table")

	);

	config.parse(argc, argv);

	// --------------------------------------------------------------------

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.has("help"))
	{
		std::cout << config << std::endl
			 << std::endl
			 << R"(Tortoize validates protein structure models by checking the 
Ramachandran plot and side-chain rotamer distributions. Quality
Z-scores are given at the residue level and at the model level 
(ramachandran-z and torsions-z). Higher scores are better. To compare 
models or to describe the reliability of the model Z-scores jackknife-
based standard deviations are also reported (ramachandran-jackknife-sd 
and torsion-jackknife-sd).

References: 
- Sobolev et al. A Global Ramachandran Score Identifies Protein 
  Structures with Unlikely Stereochemistry, Structure (2020),
  DOI: https://doi.org/10.1016/j.str.2020.08.005
- Van Beusekom et al. Homology-based loop modeling yields more complete
  crystallographic  protein structures, IUCrJ (2018),
  DOI: https://doi.org/10.1107/S2052252518010552
- Hooft et al. Objectively judging the quality of a protein structure
  from a Ramachandran plot, CABIOS (1993),
  DOI: https://doi.org/10.1093/bioinformatics/13.4.425 
)" << std::endl
			 << std::endl;
		exit(0);
	}
	
	if (config.has("build"))
	{
		buildDataFile(config.get<std::string>("build"));
		exit(0);
	}

	if (config.operands().empty())
	{
		std::cerr << "Input file not specified" << std::endl;
		exit(1);
	}

	cif::VERBOSE = config.count("verbose");

	if (config.has("log"))
	{
		if (config.operands().size() != 2)
		{
			std::cerr << "If you specify a log file, you should also specify an output file" << std::endl;
			exit(1);
		}

		std::string logFile = config.get<std::string>("log");
		
		// open the log file
		int fd = open(logFile.c_str(), O_CREAT|O_RDWR, 0644);
		if (fd < 0)
			throw std::runtime_error("Opening log file " + logFile + " failed: " + strerror(errno));
	
		// redirect stdout and stderr to the log file
		dup2(fd, STDOUT_FILENO);
		dup2(fd, STDERR_FILENO);
		close(fd);
	}

	if (config.has("dict"))
	{
		for (auto dict: config.get<std::vector<std::string>>("dict"))
			cif::compound_factory::instance().push_dictionary(dict);
	}

	// --------------------------------------------------------------------
	
	json data = tortoize_calculate(config.operands().front());

	if (config.has("output"))
	{
		std::ofstream of(config.get<std::string>("output"));
		if (not of.is_open())
		{
			std::cerr << "Could not open output file" << std::endl;
			exit(1);
		}
		of << data;
	}
	else
		std::cout << data << std::endl;
	
	return 0;
}

// --------------------------------------------------------------------

// recursively print exception whats:
void print_what (const std::exception& e)
{
	std::cerr << e.what() << std::endl;
	try
	{
		std::rethrow_if_nested(e);
	}
	catch (const std::exception& nested)
	{
		std::cerr << " >> ";
		print_what(nested);
	}
}

int main(int argc, char* argv[])
{
	int result = -1;
	
	try
	{
		result = pr_main(argc, argv);
	}
	catch (std::exception& ex)
	{
		print_what(ex);
		exit(1);
	}

	return result;
}
