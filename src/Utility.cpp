// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Lance Hepler

#include <climits>
#include <cstdlib>
#include <fstream>

#include <sys/stat.h>

#include <boost/algorithm/string.hpp>

#include <pacbio/consensus/Template.h>
#include <pacbio/io/Utility.h>

using namespace std;

namespace PacBio {
namespace IO {
namespace {

// Forward declare Recursor for WriteMatrix
class AbstractRecursor;

std::ostream& operator<<(std::ostream& out, const std::pair<size_t, size_t>& x)
{
    return out << '(' << x.first << ", " << x.second << ')';
}

}  // anonymous namespace

string AbsolutePath(const string& path)
{
    char buf[PATH_MAX + 1];
    char* res = realpath(path.c_str(), buf);

    if (res != nullptr) {
        string abspath(res);
        return abspath;
    }

    return path;
}

bool FileExists(const string& path)
{
    struct stat buffer;
    return stat(path.c_str(), &buffer) == 0;
}

string FileExtension(const string& path)
{
    size_t fileStart = path.find_last_of("/");

    if (fileStart == string::npos) fileStart = 0;

    // increment beyond the '/'
    ++fileStart;

    size_t extStart = path.substr(fileStart, path.length() - fileStart).find_last_of(".");

    if (extStart == string::npos) return "";

    // increment beyond '.'
    ++extStart;

    return path.substr(fileStart + extStart, path.length() - fileStart - extStart);
}

void FlattenFofn(vector<string>& res, const string& file)
{
    using boost::algorithm::iends_with;
    using boost::algorithm::trim;

    if (iends_with(file, ".fofn")) {
        ifstream fofn(file);
        string line;
        while (getline(fofn, line)) {
            trim(line);
            FlattenFofn(res, line);
        }
    } else if (iends_with(file, ".bam"))
        res.push_back(file);
    else
        throw invalid_argument("not a .fofn or .bam file!");
}

vector<string> FlattenFofn(const vector<string>& files)
{
    vector<string> res;

    for (const auto& file : files) {
        FlattenFofn(res, file);
    }

    return res;
}

bool ValidBaseFeatures(const PacBio::BAM::DataSet& ds)
{   
    for (const auto& bam : ds.BamFiles()) {
        for (const auto& rg : bam.Header().ReadGroups()) {
            // P6-C4 and S/P1-C1/beta do not require covariates besides SNR
            if (rg.SequencingChemistry() == "P6-C4" || rg.SequencingChemistry() == "S/P1-C1/beta")
                continue;
            // everything else requires IPD and PulseWidth
            else if (!rg.HasBaseFeature(PacBio::BAM::BaseFeature::IPD) ||
                     !rg.HasBaseFeature(PacBio::BAM::BaseFeature::PULSE_WIDTH))
                return false;
        }
    }
    return true;
}

template <typename M>
void WriteMatrix(const M& mat, const AbstractRecursor& recursor, const size_t offset, const bool reversed)
{
    std::cout << std::pair<size_t, size_t>(mat.Rows(), mat.Columns()) << std::endl;

    for (size_t j = 0; j < mat.Columns(); ++j)
        std::cout << " " << mat.UsedRowRange(j);
    std::cout << std::endl;

    std::cout << "lg: ";
    for (size_t j = 0; j < mat.Columns(); ++j)
        std::cout << "\t" << std::fixed << std::setprecision(3) << mat.GetLogScale(j);
    std::cout << std::endl;

    std::cout << "lgS: ";
    double lgS = 0.0;
    for (size_t j = 0; j < mat.Columns(); ++j)
        std::cout << "\t" << std::fixed << std::setprecision(3) << (lgS += mat.GetLogScale(j));
    std::cout << std::endl;

    for (size_t i = 0; i < mat.Rows(); ++i)
    {
        for (size_t j = 0; j < mat.Columns(); ++j)
        {
            size_t col = reversed ? mat.Columns() - j - 1 : j;
            std::cout << "\t" << std::fixed << std::setprecision(3) << std::log(mat.Get(i, j)) + mat.GetLogScale(j) + recursor.UndoCounterWeights(col + offset);
        }
        std::cout << std::endl;
    }
}

}  // namespace CCS
}  // namespace PacBio
