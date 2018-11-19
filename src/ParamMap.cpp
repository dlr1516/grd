#include <grd/ParamMap.h>
#include <cstdio>
#include <sstream>

bool ParamMap::read(std::istream& in) {
    std::string name;
    std::string value;
    std::string length;
    char buffer[MAX_LEN];
    while (in.getline(buffer, MAX_LEN)) {
        std::stringstream ss(buffer);
        if (ss >> name >> value) {
            if (name.empty() || name.at(0) == COMMENT) {
                continue;
            }
            table_.erase(name);
            table_.insert(std::make_pair(name, value));
        }
    }
    //  while (in >> name) {
    //    if (name.at(0) == COMMENT) {
    //      in.getline(buffer, MAX_LEN);
    //      continue;
    //    }
    //    if (in >> value) {
    //      table_.erase(name);
    //      table_.insert(std::make_pair(name, value));
    //    }
    //  }
    return true;
}

bool ParamMap::read(std::string& filename) {
    std::ifstream file(filename.c_str());
    if (!file) {
        std::cerr << __FILE__ << "," << __LINE__ << ": Cannot open \"" << filename << "\"" << std::endl;
        return false;
    }
    read(file);
    file.close();
    return true;
}

bool ParamMap::read(int argc, char** argv) {
    int argi = 1;
    while (argi < argc) {
        //        std::cout << "parameter " << argi << std::endl;
        std::string paramName(argv[argi]);
        if (argi + 1 < argc && isOption(paramName)) { // Checks if paramName is an option
            std::string paramNameStripped = paramName.substr(1);
            std::string paramValue(argv[argi + 1]);
            if (!isOption(paramValue)) {
                //                std::cout << "Insert variable \"" << paramNameStripped << "\", value \"" << paramValue << "\"" << std::endl;
                table_.erase(paramNameStripped);
                table_.insert(std::make_pair(paramNameStripped, paramValue));
                //                const_iterator v = table_.find(paramNameStripped);
                //                std::cout << "after find" << std::endl;
                //                if (v != table_.end()) {
                //                    std::cout << "Table: find(" << v->first << "), value \"" << v->second << "\"" << std::endl;
                //                } else {
                //                    std::cout << "Table: find(" << paramNameStripped << "), NO VALUE" << std::endl;
                //                }
            }
        }
        ++argi;
    }
}

bool ParamMap::write(std::ostream& out) const {
    return write(out, "");
}

bool ParamMap::write(std::ostream& out, std::string prefix) const {
    for (const_iterator it = table_.begin(); it != table_.end(); ++it) {
        out << prefix << it->first << " " << it->second << std::endl;
    }
    return true;
}

bool ParamMap::write(std::string& filename) const {
    return write(filename, "");
}

bool ParamMap::write(std::string& filename, std::string prefix) const {
    std::ofstream file(filename.c_str());
    if (!file) {
        std::cerr << __FILE__ << "," << __LINE__ << ": Cannot open \"" << filename << "\"" << std::endl;
        return false;
    }
    write(file, prefix);
    file.close();
    return true;
}

void ParamMap::setParamString(std::string paramName, std::string paramValue) {
    table_.erase(paramName);
    table_.insert(std::make_pair(paramName, paramValue));
}

bool ParamMap::isOption(std::string str) {
    return (str.length() >= 2 && str.at(0) == '-' && isalpha(str.at(1)));
}


