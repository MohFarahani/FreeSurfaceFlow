
#include "FileReader.hh"
#include <cstdlib>



void FileReader::registerIntParameter(const std::string &key, int init)
{
    PROG("register int parameter " + key);
    ints[key] = init; // register in the right storage
}

void FileReader::registerRealParameter(const std::string &key, real init)
{
    PROG("register real parameter " + key);
    reals[key] = init; // register in the right storage
}

void FileReader::registerStringParameter(const std::string &key, const std::string &init)
{
    PROG("register string parameter " + key);
    strings[key] = init; // register in the right storage
}

void FileReader::setParameter(const std::string &key, const std::string &in)
{
    bool found = false;
    // go trough the right storage
    using namespace std;

    for (std::map<std::string, std::string>::iterator it = strings.begin(); it != strings.end(); ++it)
    {
        if ((it->first).compare(key) == 0)    // the right key is found
        {
            it->second = in; // set the value
            found = true;
        }
    }
    if (!found)
        WARN("key " + key + " was not registered!")
    }

void FileReader::setParameter(const std::string &key, real in)
{
    bool found = false;
    // go trough the right storage
    for (std::map<std::string, real>::iterator it = reals.begin(); it != reals.end(); ++it)
    {
        if (it->first == key)    // the right key is found
        {
            it->second = in; // set the value
            found = true;
        }
    }
    if (!found)
        WARN("BLINKkey " + key + " was not registered!")
    }

void FileReader::setParameter(const std::string &key, int in)
{
    bool found = false;
    // go trough the right storage
    for (std::map<std::string, int>::iterator it = ints.begin(); it != ints.end(); ++it)
    {
        if (it->first == key)    // the right key is found
        {
            it->second = in; // set the value
            found = true;
        }
    }
    if (!found)
        WARN("key " + key + " was not registered!")
    }

// helper function: check if value is a real
bool isReal(const std::string &value)
{
    const char *str = value.c_str();
    char *endptr = 0;
    strtod(str, &endptr);

    if (*endptr != '\0' || endptr == str)
        return false;
    return true;
}

// helper function: check if value is an int
bool isInt(const std::string &value)
{
    std::size_t found = value.find_last_not_of("0123456789");
    if (found == std::string::npos)   // value is an int
        return true;

    // value is not an int
    return false;
}

bool FileReader::readFile(const std::string &name)
{
    // The file (name) has to be in the same directory as this programm!
    // "FileReaderTestInput.txt" for example was moved to this directory.
    bool res = false; // assume that the reading process failed
    std::ifstream objectFile(name.c_str()); // create an object for the file
    PROG("trying to open " + name);


    // check if file is open
    if (objectFile.is_open())
    {
        PROG("open " + name);

        // read one line after another
        for (std::string line; getline(objectFile, line);)
        {
            // get text before '#' and ignore rest of the line
            std::string pline = line.substr(0, line.find("#"));

            std::size_t found = pline.find_last_not_of(" \t\f\v\n\r");
            if (found != std::string::npos)
            {
                pline.erase(found + 1);

                // print the relevant data from line
                if (!pline.empty())
                {

                    // find key and value in line
                    std::string key = pline.substr(0, pline.find_first_of(" \t\f\v\n\r"));
                    std::size_t begin = pline.find_last_of(" \t\f\v\n\r") + 1;
                    std::string value = pline.substr(begin, pline.length() - 1);

                    if (isInt(value))     // value is an int
                    {
                        setParameter(key, atoi(value.c_str()));     // set value
                    }
                    else if (isReal(value))       // value is a real
                    {
                        setParameter(key, atof(value.c_str()));     // set value
                    }
                    else     // value is a string
                    {
                        setParameter(key, value);   // set value
                    }
                }
            }
        }
        // reading was successful
        res = true;
        // close file
        PROG("close " + name);
        objectFile.close();
    }

    return res;
}



void FileReader::printParameters() const
{
    // go trough every storage and print its input
    for (std::map<std::string, int>::const_iterator it = ints.begin(); it != ints.end(); ++it)
        std::cout << it->first << " " << it->second << std::endl;

    for (std::map<std::string, real>::const_iterator it = reals.begin(); it != reals.end(); ++it)
        std::cout << it->first << " " << it->second << std::endl;

    for (std::map<std::string, std::string>::const_iterator it = strings.begin(); it != strings.end(); ++it)
        std::cout << it->first << " " << it->second << std::endl;
}


