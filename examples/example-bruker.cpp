#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdio>
struct ProgramOptions 
{
  double mzWindowLow_;
  double snWindowLow_;
  double mzWindowHigh_;
  double snWindowHigh_;
  std::string inputfileName_;
  std::string outputfileName_;
};


struct Centroid 
{
  double rt_;
  double mz_;
  double sn_;
  double abundance_;
  Centroid(const double &rt, const double& mz, const double & sn, 
    const double & abundance) 
    : rt_(rt), mz_(mz), sn_(sn), abundance_(abundance) {}
};

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


std::vector<Centroid>
parseFile(ProgramOptions & options) {
  std::vector<Centroid> centroids;
  std::ifstream ifs(options.inputfileName_.c_str());
  std::string str;
  double mz, rt, sn, abundance;
  std::string pol, mode, mslevel, line, massrange;
  int unknown, numentries;
  typedef boost::tokenizer<boost::char_separator<char> > Tokenizer;
  boost::char_separator<char> sep(",");

  while(getline(ifs, str)) {
    Tokenizer tokens(str, sep);
    Tokenizer::iterator it = tokens.begin();
    sscanf(it->c_str(), "%e", &rt);
    ++it;
    sscanf(it->c_str(), "%s", &pol);
    std::cout << rt << " " << pol << std::endl;


  }




}

int main(int argc, char * argv[]) {

  ProgramOptions options;
  if (!parseProgramOptions(argc, argv, options)) {
    return 0;
  }
  std::vector<Centroid> centroids = parseFile(options);


}
