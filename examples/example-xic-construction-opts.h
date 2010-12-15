/* $Id: example-xic-construction-opts.h 1 2010-10-30 01:14:03Z mkirchner $
 *
 * Copyright (c) 2010 Buote Xu <buote.xu@gmail.com>
 * Copyright (c) 2010 Marc Kirchner <marc.kirchner@childrens.harvard.edu>
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

#ifndef __LIBFBI_EXAMPLES_PROGRAMPARSER_H__
#define __LIBFBI_EXAMPLES_PROGRAMPARSER_H__

#include <boost/program_options.hpp>
#include "example-xic-construction.h"

int parseProgramOptions(int argc, char* argv[], ProgramOptions& options)
{
  namespace po = boost::program_options;
  std::string config_file;
  po::options_description generic("Generic options");
  generic.add_options()
    ("help", "Display this help message")
    ("config,c", po::value<std::string>(&config_file), "config file")
    ("inputfile,i", po::value<std::string>(&options.inputfileName_), "input file")
    ("outputfile,o", po::value<std::string>(&options.outputfileName_), "output file")
    ;

  po::options_description config("Allowed options");
  config.add_options()
      ("mzLow", po::value<double>(&options.mzWindowLow_)->default_value(
        std::numeric_limits<double>::min()),
        "Lower bound of the query window for m/z")
      ("snLow", po::value<double>(&options.snWindowLow_)->default_value(
        std::numeric_limits<double>::min()),
        "Lower bound of the query window for the scan number")
      ("mzHigh", po::value<double>(&options.mzWindowHigh_)->default_value(
        std::numeric_limits<double>::max()),
        "Upper bound of the query window for m/z")
      ("snHigh", po::value<double>(&options.snWindowHigh_)->default_value(
        std::numeric_limits<double>::max()),
         "Upper bound of the query window for the scan number");

  po::options_description cmdline_options("Options available via command line");
  cmdline_options.add(generic).add(config);

  po::options_description config_file_options(
    "Options available in the config file");
  config_file_options.add(config);

  po::options_description visible("Allowed options");
  visible.add(generic).add(config);
  
  po::positional_options_description p;
  p.add("inputfile", -1);
  
  po::variables_map vm;

  po::store(po::command_line_parser(argc, argv).options(
    cmdline_options).positional(p).run(), vm);
  po::notify(vm);    
  if (vm.count("config")){ 
    std::ifstream ifs(config_file.c_str());
    if (!ifs) {
      std::cout << "can't open config file: " << config_file << '\n';
      return 0;
    } else {
      std::cerr << "Using config file" << config_file << '\n';
      po::store(po::parse_config_file(ifs, config_file_options), vm);
      po::notify(vm);    
    }
  }

  if (vm.count("help")) {
    std::cout << visible << "\n";
    return 0;
  }
  
  if (!vm.count("inputfile")) {
    std::cerr << "InputFile needed" << '\n';
    std::cout << visible << '\n';
    return 0;
  }

  if (!vm.count("outputfile")) {
    options.outputfileName_ = options.inputfileName_ + std::string(".out");
  }
  return 1;
}

#endif
