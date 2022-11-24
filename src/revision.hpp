// Generated revision file

#pragma once

#include <ostream>

const char kProjectName[] = "tortoize";
const char kVersionNumber[] = "2.0.9";
const char kVersionGitTag[] = "8507ac2";
const char kBuildInfo[] = "110*";
const char kBuildDate[] = "2022-11-24T10:26:13Z";

inline void write_version_string(std::ostream &os, bool verbose)
{
	os << kProjectName << " version " << kVersionNumber << std::endl;
	if (verbose)
	{
		os << "build: " << kBuildInfo << ' ' << kBuildDate << std::endl;
		if (kVersionGitTag[0] != 0)
			os << "git tag: " << kVersionGitTag << std::endl;
	}
}
