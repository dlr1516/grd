#include <grd/fileutils.h>
#include <glob.h>
#include <sstream>
#include <boost/date_time/posix_time/posix_time.hpp>

// ----------------------------------------------
// I/O OPERATIONS
// ----------------------------------------------

void glob(const std::string globPath, std::vector<std::string>& matchingFiles) {
    glob_t glob_result;
    matchingFiles.clear();
    ::glob(globPath.c_str(), GLOB_TILDE, NULL, &glob_result);
    for (unsigned int i = 0; i < glob_result.gl_pathc; ++i) {
        matchingFiles.push_back(std::string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
}

std::string generateStampedString(const std::string prefix, const std::string postfix) {
    boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
    std::ostringstream formatter;
    std::string formatstring = prefix + "%Y%m%d_%H%M_%S" + postfix;
    formatter.imbue(std::locale(std::cout.getloc(), new boost::posix_time::time_facet(formatstring.c_str())));
    formatter << now;
    return formatter.str();
}

std::string getWithoutExtension(std::string filename) {
    boost::filesystem::path path(filename);
    return path.stem().string();
}

std::string getPrefix(std::string filename) {
    // Strips filename of the path 
    boost::filesystem::path filepath(filename);
    std::string name = filepath.filename().string();
    std::string prefix;
    //  std::cout << "  name: \"" << name << "\"\n";

    // Finds the prefix
    size_t pos = name.find_first_of('_');
    if (pos != std::string::npos) {
        prefix = name.substr(0, pos);
    } else {
        prefix = name;
    }
    return prefix;
}

std::string getShortName(std::string filename) {
    std::stringstream ss;
    std::string prefix = getPrefix(filename);
    boost::filesystem::path filenamePath = filename;
    filename = filenamePath.filename().string();
    // Computes a digest on the string
    unsigned int h = 19;
    for (int i = 0; i < filename.length(); ++i) {
        h = ((h * 31) + (unsigned int) filename[i]) % 97;
    }
    //  std::cout << "\nglob \"" << filenamePath.string() << "\" filename \"" << filename << "\" hash " << h << std::endl;
    ss << prefix << "_" << std::setw(2) << std::setfill('0') << h;
    return ss.str();
}

std::string getLeafDirectory(std::string filename) {
    boost::filesystem::path filenamePath = filename;
    std::string parent = filenamePath.parent_path().string();
    size_t pos = parent.find_last_of('/');
    std::string leafDir = "";
    if (pos != std::string::npos) {
        leafDir = parent.substr(pos + 1, parent.length());
    }
    return leafDir;
}

