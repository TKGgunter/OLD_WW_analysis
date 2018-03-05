//Thoth Gunter
#include <bitset>
#include <iostream>
#include <stdint.h>

int main()
{
    float f = 0.125;
    char* bits = reinterpret_cast<char*>(&f);
    for(std::size_t n = 0; n < sizeof f; ++n)
            std::cout << std::bitset<8>(bits[n]);
    std::cout << '\n';
    std::cout << bits << '\n';



    int i; // declare variables
    f = 0.125;
    int bit = 0;
    std::cout << "Enter a floating point number: "; // enter a float

    int *b = reinterpret_cast<int*>(&f); // use reinterpret_cast function
    for (int k = 31; k >=0; k--) // for loop to print out binary pattern
    {
      bit = ((*b >> k)&1); // get the copied bit value shift right k times, then and with a 1.
      std::cout << bit; // print the bit.
    }


  std::cout << std::endl; 
  std::string str ("2\x00\x16\x00'\x00\x02\x00)\x00\x03\x00-\x00\x04\x00\x00\x00\x15\x00F\x00", 22);
  for (int i=0; i<str.length()/2; ++i)
  {
    std::cout << ((uint16_t)str[i*2] | (str[i*2+1] << 8)) << " from " << (uint16_t) str[i*2] << " and " << (uint16_t)  str[i*2 + 1] <<"\n";
  }
  std::cout << std::endl; 



}


