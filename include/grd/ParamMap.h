/**
 * GRD - Geometric Relation Distribution
 * Copyright (C) 2018 Dario Lodi Rizzini.
 *
 * GRD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * GRD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with GRD.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef GRD_PARAM_MAP_H
#define GRD_PARAM_MAP_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <boost/lexical_cast.hpp>

namespace grd {

    /** Reads and stores parameters from a string, a file, etc.
     */
    class ParamMap {
    public:
        //typedef std::unordered_map<std::string, std::string> table_type;
        typedef std::map<std::string, std::string> table_type; // a map stores lexically ordered parameters (nicer to view!)
        typedef table_type::iterator iterator;
        typedef table_type::const_iterator const_iterator;

        const static char COMMENT = '#';
        const static unsigned int MAX_LEN = 2000;

        /** Default constructor.
         */
        ParamMap() : table_() {
        }

        /** Destructor.
         */
        ~ParamMap() {
        }

        /** Clears all the content of param table.
         */
        void clear() {
            table_.clear();
        }

        /** Reads params from an input stream in the format:
         *   key1 value1
         *   key2 value2
         *   ...
         */
        bool read(std::istream& in);

        /** Reads params from an input file (format as above).
         */
        bool read(std::string& filename);

        /** Reads from a command line. Required format
         *   ...
         *   argv[i] = "-key1"      // starting with "-"
         *   argv[i+1] = "value1"
         */
        bool read(int argc, char** argv);

        /** Writes the pairs (key,value).
         */
        bool write(std::ostream& out) const;

        /**
         * Writes the pairs (key,value) by displaying a prefix before each pair. 
         * @param out
         * @param prefix
         * @return 
         */
        bool write(std::ostream& out, std::string prefix) const;

        /** Writes the parameters to an output file (format as above).
         */
        bool write(std::string& filename) const;

        /**
         * Writes the parameters to an output file (format as above) by displaying 
         * a prefix before each parameter pair. 
         * @param filename
         * @param prefix
         * @return 
         */
        bool write(std::string& filename, std::string prefix) const;

        /** Sets the param (as a set of string).
         */
        void setParamString(std::string paramName, std::string paramValue);

        /** Sets the param (as a set of string).
         */
        template <typename Value>
        void setParam(std::string paramName, const Value& paramValue) {
            std::stringstream sstr;
            sstr << paramValue;
            table_.erase(paramName);
            table_.insert(std::make_pair(paramName, sstr.str()));
        }

        /** Casts the string value of a given parameters to the desired value.
         */
        template <typename Value>
        bool getParam(std::string paramName, Value& value, const Value& defaultValue) {
            const_iterator v = table_.find(paramName);
            if (v != table_.end()) {
                try {
                    value = boost::lexical_cast<Value>(v->second);
                } catch (boost::bad_lexical_cast const&) {
                    std::cerr << __FILE__ << "," << __LINE__ << ": Error: cannot cast string \""
                            << v->second << "\" to type \"" << typeid (Value).name() << "\" for variable \"" << v->first << "\"" << std::endl;
                }
            } else {
                //            std::cerr << "Parameter " << paramName << " not found." << std::endl;
                value = defaultValue;
                setParam(paramName, defaultValue);
                return false;
            }
        }

    protected:
        table_type table_;

        static bool isOption(std::string str);
    };

} // end of namespace 

#endif