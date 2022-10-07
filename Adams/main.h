#ifndef MAIN_H
#define MAIN_H
inline const char* VERSION = "Version:\n  2.2 (2022-10-06)";

//#define MYDEPLOY

int main_res(int argc, char** argv, int index);
int main_prod(int argc, char** argv, int index);
int main_export(int argc, char** argv, int index);

int main_2cell(int argc, char** argv, int index);
int main_2cell_export(int argc, char** argv, int index);

#endif