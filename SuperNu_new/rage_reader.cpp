/*
 * main.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: dmarce1
 */

#include <stdio.h>
#include <string>
#include <cstring>
#include <vector>
#include <list>
#include <map>
#include <algorithm>

#define nrage_fracs 17
const char* rage_fracs[nrage_fracs] = { "h", "he", "c", "n", "o", "ne", "mg", "si", "s", "ar", "ca", "ti", "cr", "fe", "ni", "ni56",
		"fe52" };

std::string read_line(FILE* fp) {
	std::string iline;
	if (feof(fp)) {

	} else {
		char c = fgetc(fp);
		while (c != '\n' && !feof(fp)) {
			iline += c;
			c = fgetc(fp);
		}
	}
	return iline;
}

std::vector<std::string> get_tokens(const std::string& iline) {
	const char* ptr = iline.c_str();
	std::vector<std::string> tokens;
	while (*ptr != '\0' && *ptr != '\n') {
		std::string token;
		while ((isspace(*ptr) || *ptr == '#') && *ptr != '\0') {
			++ptr;
		}
		while (!isspace(*ptr) && *ptr != '\0') {
			token += *ptr;
			++ptr;
		}
		tokens.push_back(token);
	}
	return tokens;
}

void read_file(const char* filename, std::map<std::string, std::vector<double> >& data_map,
		std::vector<std::string>& field_ids) {
	FILE* fp = fopen(filename, "rt");
	int ncols = 0;
	int nrows = 0;
	int line_cnt = 0;
	if (fp != NULL) {
		std::vector<std::string> tokens;
		auto iline = read_line(fp);
		if (iline.size()) {
			tokens = get_tokens(iline);
			ncols = tokens.size();
		}
		if (ncols != 0) {
			fprintf(stderr, "Reading %i columns of data\n", ncols);
			fprintf(stderr, "Field identifiers: \n");
			for (int i = 0; i != ncols; ++i) {
				field_ids.push_back(tokens[i]);
				fprintf(stderr, "%s\n", tokens[i].c_str());
			}
			while (!feof(fp)) {
				auto iline = read_line(fp);
				line_cnt++;
				if (iline.size()) {
					tokens = get_tokens(iline);
					if (int(tokens.size()) != ncols) {
						printf("Read %i but expected %i on line %i", int(tokens.size()), ncols, line_cnt);
						break;
					} else {
						for (int i = 0; i != ncols; ++i) {
							double value = atof(tokens[i].c_str());
							data_map[field_ids[i]].push_back(value);
						}
					}
					++nrows;
				} else {
					break;
				}
			}
			for (int i = 0; i != nrage_fracs; ++i) {
				if (data_map.find(rage_fracs[i]) == data_map.end()) {
					data_map[rage_fracs[i]] = std::vector<double>(nrows, 0.0);
					field_ids.push_back(rage_fracs[i]);
				}
			}
		} else {
			printf("Empty file\n");
		}
		fclose(fp);
	} else {
		printf("%s not found\n", filename);
	}

}

const int niso = 6;
const std::pair<const char*, const char*> isopairs[6] = { std::make_pair<const char*, const char*>("ni56", "ni"),
		std::make_pair<const char*, const char*>("co56", "co"), std::make_pair<const char*, const char*>("fe52", "fe"),
		std::make_pair<const char*, const char*>("mn52", "mn"), std::make_pair<const char*, const char*>("cr48", "cr"),
		std::make_pair<const char*, const char*>("v48", "v48") };

const int nsupernu_fracs = 31;

void convert_file(std::map<std::string, std::vector<double>>& map, int N) {
	const std::vector<double>& T = map["temperature"];
	const std::vector<double>& rho = map["density"];
	const std::vector<double>& r = map["radius"];
	const std::vector<double>& vel = map["velocity"];
	double vmax, rmax;
	vmax = *(std::max_element(vel.begin(), vel.end()));
	rmax = *(std::max_element(r.begin(), r.end()));
	fprintf(stderr, "Maximum Velocity = %e\n", vmax);
	fprintf(stderr, "Maximum Radius = %e\n", rmax);
	const int nrows = r.size();

	const double dvel = vmax / (double(N));
	double t = 0.0;

	for (int f = 0; f != nrage_fracs; ++f) {
		for (int i = 0; i != niso; ++i) {
			const auto key2 = rage_fracs[f];
			if (std::strcmp(key2, isopairs[i].second) == 0) {
				const auto key1 = isopairs[i].first;
				if (map.find(key1) != map.end()) {
					for (int j = 0; j != nrows; ++j) {
						map[key2][j] += map[key1][j];
					}
				}
			}
		}
	}

	std::vector<double> mass(N, 0.0);
	std::vector<double> temp(N, 0.0);
	std::vector<std::vector<double>> frac_mass(nrage_fracs, std::vector<double>(N, 0.0));
	double mtot = 0.0;
	for (int i = 0; i != nrows; ++i) {
//		printf("%e %e %e\n", r[i], vel[i], rho[i]);
		const int vbin = vel[i] / dvel;
		const double dr0 = i == 0 ? r[i] : r[i] - r[i - 1];
		const double dr1 = i == nrows - 1 ? dr0 : r[i + 1] - r[i];
		const double dr = 0.5 * (dr0 + dr1);
		const double dV = 4.0 * M_PI * r[i] * r[i] * dr;
		const double dm = rho[i] * dV;
		mtot += dm;
		mass[vbin] += dm;
		temp[vbin] += T[i] * dm;
		for (int f = 0; f != nrage_fracs; ++f) {
			frac_mass[f][vbin] += map[rage_fracs[f]][i] * dm;
		}
		t += r[i] / vel[i] * dm;
	}
	t /= mtot;
	for (int vbin = 0; vbin != N; ++vbin) {
		temp[vbin] /= mass[vbin];
		for (int f = 0; f != nrage_fracs; ++f) {
			frac_mass[f][vbin] /= mass[vbin];
		}
	}
	const int nfracs = frac_mass.size();
	printf("# spherical\n");
	printf("# %i 1 1 %i %i\n", N, nfracs + 3, nfracs);
	printf("#");
	printf("%10s", "rightvel");
	printf("%12s", "mass");
	printf("%12s", "temp");
	for (int f = 0; f != nrage_fracs; ++f) {
		printf("%12s", rage_fracs[f]);
	}
	printf("\n");
	fprintf(stderr, "Time at %e s  (%e d)\n", t, t / (3600.0 * 24.0));
	for (int i = 0; i != N; ++i) {
		const double vel = dvel * double(i + 1);
		printf(" %.4e  %.4e  %.4e", vel, mass[i], temp[i]);
//		double n = 0.0;
		for (int f = 0; f != nrage_fracs; ++f) {
	//		n += frac_mass[f][i];
			printf("  %.4e", frac_mass[f][i]);
		}
	//	fprintf(stderr, "%e\n", 1.0 - n);
		printf("\n");
	}
}

int main(int argc, char* argv[]) {
	std::map<std::string, std::vector<double> > data_map;
	std::vector<std::string> field_ids;
	if (argc == 3) {
		read_file(argv[1], data_map, field_ids);
		int N = atoi(argv[2]);
		convert_file(data_map, N);
	} else {
		printf("arguments <filename> <nrows> \n");
	}
	return 0;
}
