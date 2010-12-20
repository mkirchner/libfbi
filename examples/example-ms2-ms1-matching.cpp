/* 
 * example-ms2-ms1-matching.cpp
 *
 * Copyright (c) 2010 Marc Kirchner <marc.kirchner@childrens.harvard.edu>
 * Copyright (c) 2010 Buote Xu <buote.xu@gmail.com>
 *
 * This file is part of libfbi.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without  restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions: 
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR  OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/time.h>
#include <tuple>
#include <utility>
#include <vector>

#include "fbi/tuplegenerator.h"
#include "fbi/fbi.h"
#include "fbi/connectedcomponents.h"

/*
 * User classes
 */
struct Xic
{
    double rt_;
    double mz_;
    double abundance_;
    Xic(const double &rt, const double& mz, const double & abundance) 
      : rt_(rt), mz_(mz), abundance_(abundance) {}
};

struct MS2Scan
{
    double rt_;
    double mz_;
    typedef std::vector<std::pair<double, double> > Ions;
    Ions ions_;
    MS2Scan(const double &rt, const double& mz, const Ions& ions)
      : rt_(rt), mz_(mz), ions_(ions) {}
};

/*
 * Type traits info for libfbi
 */
namespace fbi {

template<> struct Traits<Xic> : mpl::TraitsGenerator<double, double> {};
template<> struct Traits<MS2Scan> : mpl::TraitsGenerator<double, double> {};

} 

/*
 * Xic adapters
 */
struct XicBoxGenerator
{
  template <size_t N>
  typename std::tuple_element<N, 
    typename fbi::Traits<Xic>::key_type>::type 
  get(const Xic &) const;

  double fullScanPpm_;
  double rtWindow_;

  XicBoxGenerator(double fullScanPpm, double rtWindow)
    : fullScanPpm_(fullScanPpm), rtWindow_(rtWindow)
  {}
};

template <>
std::pair<double, double>  
XicBoxGenerator::get<0>(const Xic& xic) const 
{
  return std::make_pair(
    xic.mz_* (1 - fullScanPpm_ * 1E-6), 
    xic.mz_* (1 + fullScanPpm_ * 1E-6));
}

template <>
std::pair<double, double>  
XicBoxGenerator::get<1>(const Xic & xic) const
{
  return std::make_pair(xic.rt_ - rtWindow_, xic.rt_ + rtWindow_);
}

/*
 * MS2Scan classes and adapters
 */
struct MS2ScanBoxGenerator
{
  template <size_t N>
  typename std::tuple_element<N, 
    typename fbi::Traits<MS2Scan>::key_type>::type 
  get(const MS2Scan &) const;

  double preScanPpm_;
  double rtWindow_;

  MS2ScanBoxGenerator(double preScanPpm, double rtWindow)
    : preScanPpm_(preScanPpm), rtWindow_(rtWindow)
  {}
};

template <>
std::pair<double, double>  
MS2ScanBoxGenerator::get<0>(const MS2Scan& ms2scan) const 
{
  return std::make_pair(
    ms2scan.mz_* (1 - preScanPpm_ * 1E-6), 
    ms2scan.mz_* (1 + preScanPpm_ * 1E-6));
}

template <>
std::pair<double, double>  
MS2ScanBoxGenerator::get<1>(const MS2Scan & ms2scan) const
{
  return std::make_pair(ms2scan.rt_ - rtWindow_, ms2scan.rt_ + rtWindow_);
}

/*
 * store and process user options
 */

struct ProgramOptions 
{
    double prescanPpm_;
    double fullscanPpm_;
    double rtWindow_;
    std::string xicFileName_;
    std::string ms2scanFileName_;
    std::string outputFileName_;
};

int parseProgramOptions(int argc, char* argv[], ProgramOptions& options)
{
  namespace po = boost::program_options;
  std::string config_file;
  double pres, fres;
  po::options_description generic("Generic options");
  generic.add_options()
    ("help", "Display this help message")
    ("config,c", po::value<std::string>(&config_file), "config file")
    ("xicfile,x", po::value<std::string>(&options.xicFileName_), "input file")
    ("ms2scanfile,m", po::value<std::string>(&options.ms2scanFileName_), "input file")
    ("outputfile,o", po::value<std::string>(&options.outputFileName_), "output file")
    ;

  po::options_description config("Allowed options");
  config.add_options()
      ("pr", po::value<double>(&pres)->default_value(7500.0),
        "The resolution of the prescan.")
      ("fr", po::value<double>(&fres)->default_value(60000.0),
        "The resolution of the full MS1 parent scan.")
      ("rt", po::value<double>(&options.rtWindow_)->default_value(60.0),
        "The size of the retention time search window (in s).");

  po::options_description cmdline_options("Options available via command line");
  cmdline_options.add(generic).add(config);

  po::options_description config_file_options(
    "Options available in the config file");
  config_file_options.add(config);

  po::options_description visible("Allowed options");
  visible.add(generic).add(config);
  
  po::positional_options_description p;
  p.add("xicfile", 1);
  p.add("ms2scanfile", 1);
  
  po::variables_map vm;

  po::store(po::command_line_parser(argc, argv).options(
    cmdline_options).positional(p).run(), vm);
  po::notify(vm);    
  if (vm.count("config")){ 
    std::ifstream ifs(config_file.c_str());
    if (!ifs) {
      std::cout << "can't open config file: " << config_file << '\n';
      return -1;
    } else {
      std::cerr << "Using config file" << config_file << '\n';
      po::store(po::parse_config_file(ifs, config_file_options), vm);
      po::notify(vm);    
    }
  }

  if (vm.count("help")) {
    std::cout << visible << "\n";
    return -1;
  }
  
  if (vm.count("xicfile")) {
    std::cerr << "xicfile: " <<
        options.xicFileName_ << '\n';
  } else {
    std::cerr << "XIC file required." << '\n';
    std::cout << visible << '\n';
    return -1;
  }

  if (vm.count("ms2scanfile")) {
    std::cerr << "ms2scanfile: " <<
        options.ms2scanFileName_ << '\n';
  } else {
    std::cerr << "MS2 input file required." << '\n';
    std::cout << visible << '\n';
    return -1;
  }

  if (vm.count("outputfile")) {
    std::cerr << "outputfile: " <<
        options.outputFileName_ << '\n';
  } else {
    options.outputFileName_ = options.xicFileName_ + std::string(".out");
    std::cerr << "outputfile: " <<
        options.outputFileName_ << '\n';
  }

  options.prescanPpm_ = 1e6 / (pres * (2*std::sqrt(2*std::log(2))));
  options.fullscanPpm_ = 1e6 / (fres * (2*std::sqrt(2*std::log(2))));

  return 0;
}

/*
 * load xics from file
 */
std::vector<Xic> parseXicFile(ProgramOptions& options)
{
  std::vector<Xic> xics;
  std::ifstream ifs(options.xicFileName_.c_str());
  ifs.setf(std::ios::fixed, std::ios::floatfield);

  double mz, rt, abundance;

  while (ifs >> rt >> mz >> abundance) {
    xics.push_back(Xic(rt, mz, abundance));
  }
  return xics;
}

/*
 * load MS2 scans from file
 */
std::vector<MS2Scan> parseMS2ScanFile(ProgramOptions & options)
{
  std::vector<MS2Scan> ms2scans;
  std::ifstream ifs(options.ms2scanFileName_.c_str());
  ifs.setf(std::ios::fixed, std::ios::floatfield);

  double mz, rt;
  std::string seq;

  while (ifs >> rt >> mz) {
    MS2Scan::Ions ions;
    ms2scans.push_back(MS2Scan(rt, mz, ions));
  }
  return ms2scans;
}


int main(int argc, char* argv[])
{
    using namespace fbi;

    ProgramOptions options;
    if (parseProgramOptions(argc, argv, options) != 0) {
        return -1;
    }
    // load data
    std::vector<Xic> xics = parseXicFile(options);
    std::vector<MS2Scan> ms2scans = parseMS2ScanFile(options);

    timeval start, end; 
    gettimeofday(&start, NULL);
    auto adjList = SetA<Xic, 0, 1>::SetB<MS2Scan, 0, 1>::intersect(
      xics, XicBoxGenerator(options.fullscanPpm_, options.rtWindow_),
      ms2scans, MS2ScanBoxGenerator(options.prescanPpm_, options.rtWindow_));
    gettimeofday(&end, NULL);
    std::cout << "fbi runtime: "
      << static_cast<double>(end.tv_sec - start.tv_sec) +
         static_cast<double>(end.tv_usec - start.tv_usec)* 1E-6 << '\n';

    std::ofstream ofs(options.outputFileName_.c_str());
    ofs.setf(std::ios::fixed, std::ios::floatfield);
    for (size_t i = 0; i < xics.size(); ++i) {
        Xic& xic = xics[i];
        ofs << xic.rt_ << '\t' << xic.mz_ << '\t' << '\t' 
          << xic.abundance_;
        if (adjList[i].empty()) {
            ofs << '\n';
        } else {
            ofs << " -> ";
            typedef std::set<unsigned int>::const_iterator SI;
            for (SI j = adjList[i].begin(); j != adjList[i].end(); ++j) {
                ofs << " (" << ms2scans[(*j)-xics.size()].rt_
                  << ',' << ms2scans[(*j)-xics.size()].mz_ << ')';
            }
            ofs << '\n';
       }
    }
    return 0;
}

