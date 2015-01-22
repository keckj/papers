
#ifndef PLOTUTILS_H
#define PLOTUTILS_H

namespace PlotUtils {

    template <typename T, unsigned int N>
    std::string prettyValue(T value) {
       
        std::stringstream ss;
        ss << std::fixed << std::setprecision(N) << value;

        std::string s = ss.str();

        if(N > 0 && s[s.find_last_not_of('0')] == '.') {
            s.erase(s.size() - N + 1);
        }
        return s;
    }
}


#endif /* end of include guard: PLOTUTILS_H */
