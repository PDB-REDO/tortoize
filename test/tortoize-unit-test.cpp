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

#include "tortoize.hpp"

#include <catch2/catch_all.hpp>
#include <filesystem>

namespace fs = std::filesystem;

// --------------------------------------------------------------------

// --------------------------------------------------------------------

fs::path gTestDir = fs::current_path();

int main(int argc, char *argv[])
{
	Catch::Session session; // There must be exactly one instance

	// Build a new parser on top of Catch2's
#if CATCH22
	using namespace Catch::clara;
#else
	// Build a new parser on top of Catch2's
	using namespace Catch::Clara;
#endif

	auto cli = session.cli()                               // Get Catch2's command line parser
	           | Opt(gTestDir, "data-dir")                 // bind variable to a new option, with a hint string
	                 ["-D"]["--data-dir"]                  // the option names it will respond to
	           ("The directory containing the data files") // description string for the help output
	           | Opt(cif::VERBOSE, "verbose")["-v"]["--cif-verbose"]("Flag for cif::VERBOSE");

	// Now pass the new composite back to Catch2 so it uses that
	session.cli(cli);

	// Let Catch2 (using Clara) parse the command line
	int returnCode = session.applyCommandLine(argc, argv);
	if (returnCode != 0) // Indicates a command line error
		return returnCode;

	cif::add_data_directory(gTestDir / ".." / "rsrc");

	return session.run();
}

// --------------------------------------------------------------------

TEST_CASE("first_test")
{
	auto a = tortoize_calculate(gTestDir / "1cbs.cif.gz");

	std::ifstream bf(gTestDir / "1cbs.json");

	auto b = zeep::el::object::parse_JSON(bf);

	auto ma = a["model"]["1"];
	auto mb = b["model"]["1"];

	CHECK_THAT(ma["ramachandran-jackknife-sd"].template get<double>(), Catch::Matchers::WithinRel(mb["ramachandran-jackknife-sd"].template get<double>(), 0.1));
	CHECK_THAT(ma["ramachandran-z"].template get<double>(), Catch::Matchers::WithinRel(mb["ramachandran-z"].template get<double>(), 0.1));

	CHECK_THAT(ma["torsion-jackknife-sd"].template get<double>(), Catch::Matchers::WithinRel(mb["torsion-jackknife-sd"].template get<double>(), 0.1));
	CHECK_THAT(ma["torsion-z"].template get<double>(), Catch::Matchers::WithinRel(mb["torsion-z"].template get<double>(), 0.1));

	// Residue level checks. The identifiers must match exactly, in particular
	// pdb.seqNum which is parsed from the PDB sequence number string.
	const zeep::el::object &ra = ma.at("residues");
	const zeep::el::object &rb = mb.at("residues");

	REQUIRE(ra.size() == rb.size());

	for (std::size_t i = 0; i < rb.size(); ++i)
	{
		const zeep::el::object &x = ra.at(i);
		const zeep::el::object &y = rb.at(i);

		CHECK(x.at("compID").get<std::string>() == y.at("compID").get<std::string>());
		CHECK(x.at("asymID").get<std::string>() == y.at("asymID").get<std::string>());
		CHECK(x.at("seqID").get<int>() == y.at("seqID").get<int>());
		CHECK(x.at("pdb").at("seqNum").get<int>() == y.at("pdb").at("seqNum").get<int>());
		CHECK(x.at("pdb").at("insCode").get<std::string>() == y.at("pdb").at("insCode").get<std::string>());

		// Z-scores can differ slightly between platforms, so compare them with a
		// relative tolerance and fall back to an absolute one near zero.
		const zeep::el::object &xa = x.at("ramachandran");
		const zeep::el::object &ya = y.at("ramachandran");

		CHECK(xa.at("ss-type").get<std::string>() == ya.at("ss-type").get<std::string>());
		CHECK_THAT(xa.at("z-score").get<double>(),
			Catch::Matchers::WithinRel(ya.at("z-score").get<double>(), 0.05) or
				Catch::Matchers::WithinAbs(ya.at("z-score").get<double>(), 0.001));

		// The torsion entry is only present for residues that have a chi angle.
		const zeep::el::object &xt = x["torsion"];
		const zeep::el::object &yt = y["torsion"];

		REQUIRE(xt.is_null() == yt.is_null());

		if (not yt.is_null())
		{
			CHECK(xt.at("ss-type").get<std::string>() == yt.at("ss-type").get<std::string>());
			CHECK_THAT(xt.at("z-score").get<double>(),
				Catch::Matchers::WithinRel(yt.at("z-score").get<double>(), 0.05) or
					Catch::Matchers::WithinAbs(yt.at("z-score").get<double>(), 0.001));
		}
	}
}
