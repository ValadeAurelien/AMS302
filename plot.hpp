#include <vector>
#include <string>
#include <iostream>

template<typename Tx, typename Ty>
int plot(const std::vector<Tx>& X, const std::vector<Ty>& Y, const std::string& opt = std::string(), const std::string& title = std::string())
{
    if (X.size() != Y.size()) {
	std::cout << "sizes don't match !" << std::endl;
	return -1;
    }
    FILE *pipe = popen("gnuplot -persistent", "w");
    fprintf(pipe, ("set title '" + title + "' ; plot '-' " + opt + "\n").c_str());
    for (int i=0; i<X.size(); i++)
	fprintf(pipe, "%f %f\n", (float) X.at(i), (float) Y.at(i));
    fflush(pipe);
    fclose(pipe);
    return 0;
}

template<typename Ty>
int plot(const std::vector<Ty>& Y, const std::string& opt = std::string(), const std::string& title = std::string())
{
    FILE *pipe = popen("gnuplot -persistent", "w");
    fprintf(pipe, ("set title '" + title + "' ; plot '-' " + opt + "\n").c_str());
    for (int i=0; i<Y.size(); i++)
	fprintf(pipe, "%f\n", (float) Y.at(i));
    fflush(pipe);
    fclose(pipe);
    return 0;
}
