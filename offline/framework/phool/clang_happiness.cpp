#include <iostream>

int main()
{
  std::cout << "jenkins searches our code base for .cpp files which" << std::endl;
  std::cout << "so far we do not have and clang-tidy throws a file not found error" << std::endl;
  std::cout << "But rather than removing this extension from jenkins and risking" << std::endl;
  std::cout << "that files with this extension get checked in and then are exlcuded" << std::endl;
  std::cout << "I created this dummy file which is not compiled but gives" << std::endl;
  std::cout << "clang-tidy something to sink its teeth in" << std::endl;
  return 0;
}
