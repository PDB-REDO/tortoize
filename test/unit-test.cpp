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

#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include <filesystem>

#include <cif++.hpp>

#include <zeep/json/parser.hpp>

#include "tortoize.hpp"

namespace fs = std::filesystem;

using json = zeep::json::element;

// --------------------------------------------------------------------

fs::path gTestDir = fs::current_path();

bool init_unit_test()
{
	cif::VERBOSE = 1;

	// not a test, just initialize test dir
	if (boost::unit_test::framework::master_test_suite().argc == 2)
	{
		gTestDir = boost::unit_test::framework::master_test_suite().argv[1];

		cif::add_data_directory(gTestDir / ".." / "rsrc");
	}

	// // do this now, avoids the need for installing
	// cif::addFileResource("mmcif_pdbx_v50.dic", gTestDir / ".." / "rsrc" / "mmcif_pdbx_v50.dic");

	// // initialize CCD location
	// cif::addFileResource("components.cif", gTestDir / ".." / "data" / "ccd-subset.cif");

	return true;
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(first_test)
{
	auto a = tortoize_calculate(gTestDir / "1cbs.cif.gz");

	std::ifstream bf(gTestDir / "1cbs.json");

	json b;
	zeep::json::parse_json(bf, b);

	a["software"]["version"] = "test";
	b["software"]["version"] = "test";

	// somehow a == b does not work correctly?
	// BOOST_CHECK_EQUAL(a, b);

	std::stringstream sa, sb;
	sa << a;
	sb << b;

	BOOST_CHECK_EQUAL(sa.str(), sb.str());

	if (sa.str() != sb.str())
	{
		std::ofstream out1("/tmp/a.json");
		out1 << sa.str();
		out1.close();

		std::ofstream out2("/tmp/b.json");
		out2 << sb.str();
		out2.close();
	}
}