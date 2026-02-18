/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "revision.hpp"
#include "tortoize.hpp"

#include <cif++/cif++.hpp>
#include <mcfp/mcfp.hpp>

#include <fstream>

namespace fs = std::filesystem;

using json = nlohmann::json;

#ifdef _MSC_VER
// MSVC stdlib.h definitions
#define STDIN_FILENO 0
#define STDOUT_FILENO 1
#define STDERR_FILENO 2
#else
#include <fcntl.h>
#endif

// --------------------------------------------------------------------

int pr_main(int argc, char *argv[]) {
  using namespace std::literals;

  auto &config = mcfp::config::instance();

  config.init(
      "tortoize [options] input [output]",
      mcfp::make_option("help,h", "Display help message"),
      mcfp::make_option("version", "Print version"),

      mcfp::make_option("verbose,v", "verbose output"),

      mcfp::make_option<std::string>("log", "Write log to this file"),

      mcfp::make_option<std::vector<std::string>>(
          "dict", "Dictionary file containing restraints for residues in this "
                  "specific target, can be specified multiple times."),

      mcfp::make_hidden_option<std::string>("build",
                                            "Build a binary data table")

  );

  config.parse(argc, argv);

  // --------------------------------------------------------------------

  if (config.has("version")) {
    write_version_string(std::cout, config.has("verbose"));
    exit(0);
  }

  if (config.has("help")) {
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

  if (config.has("build")) {
    buildDataFile(config.get<std::string>("build"));
    exit(0);
  }

  if (config.operands().empty()) {
    std::cerr << "Input file not specified" << std::endl;
    exit(1);
  }

  cif::VERBOSE = config.count("verbose");

  if (config.has("log")) {
    if (config.operands().size() != 2) {
      std::cerr
          << "If you specify a log file, you should also specify an output file"
          << std::endl;
      exit(1);
    }

    std::string logFile = config.get<std::string>("log");

    // open the log file
    int fd = open(logFile.c_str(), O_CREAT | O_RDWR, 0644);
    if (fd < 0)
      throw std::runtime_error("Opening log file " + logFile +
                               " failed: " + strerror(errno));

    // redirect stdout and stderr to the log file
    dup2(fd, STDOUT_FILENO);
    dup2(fd, STDERR_FILENO);
    close(fd);
  }

  if (config.has("dict")) {
    for (auto dict : config.get<std::vector<std::string>>("dict"))
      cif::compound_factory::instance().push_dictionary(dict);
  }

  // --------------------------------------------------------------------

  json data = tortoize_calculate(config.operands().front());

  if (config.operands().size() == 2) {
    std::ofstream of(config.operands().back());
    if (not of.is_open()) {
      std::cerr << "Could not open output file" << std::endl;
      exit(1);
    }
    of << data;
  } else
    std::cout << data << std::endl;

  return 0;
}

// --------------------------------------------------------------------

// recursively print exception whats:
void print_what(const std::exception &e) {
  std::cerr << e.what() << std::endl;
  try {
    std::rethrow_if_nested(e);
  } catch (const std::exception &nested) {
    std::cerr << " >> ";
    print_what(nested);
  }
}

int main(int argc, char *argv[]) {
  int result = -1;

  try {
    result = pr_main(argc, argv);
  } catch (std::exception &ex) {
    print_what(ex);
    exit(1);
  }

  return result;
}
