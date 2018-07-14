#include <vector>
#include <string>
#include <sstream>
#include <functional>
#include <algorithm>

/**
 *  \brief Generates a hash name.
 *
 *  Concatenates and hashes a vector of strings.
 */
class NativeHashGenerator {

public:

    static std::string generate(std::vector<std::string> arguments,
                                size_t world_pid = 0) {

        std::string hash_input = "";

        for (const std::string &arg: arguments ) {
            hash_input += arg;
        }

        hash_input += "_" + std::to_string(world_pid);

        std::hash<std::string> hashFunction;
        size_t hash_value = hashFunction(hash_input);

        std::ostringstream hash_str;
        hash_str << std::hex << hash_value;

        reverse(hash_input.begin(), hash_input.end());
        hash_value = hashFunction(hash_input);
        hash_str << hash_value;

        return hash_str.str();
    }

private:

    NativeHashGenerator() {}
    ~NativeHashGenerator() {}

};