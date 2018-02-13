#include <iostream>
#include <stdio.h>
#include <cassert>
#include "hashFunc.h"
#include "bloom.h"
#include <vector>
using namespace std;

int main()
{
    /*
     *    Create two hash functions
     */
    HashFunA *funa = new HashFunA();
    HashFunB * funb = new HashFunB();
    vector<HashFun*> hashfunclist;
    hashfunclist.push_back(funa);
    hashfunclist.push_back(funb);

    /*
     * Create Bloom object with two parameters :
     * size of the store array and list of hash functions
     */
    Bloom bloom(10000,hashfunclist);

    ///Add some words to bloom filter
    bloom.add("hello");
    bloom.add("world");
    bloom.add("ipad");
    bloom.add("iphone4");
    bloom.add("ipod");
    bloom.add("apple");
    bloom.add("banana");
    bloom.add("hello");

    /*
     * Test
     */
    char word[20];
    while(true)
    {
        cout<<"Please input a word : "<<endl;
        cin>>word;
        if(bloom.check(word))
        {
            cout<<"Word :"<<word<<" has been set in bloom filter."<<endl;
        }
        else
        {
            cout<<"Word :"<<word<<" not exist !" <<endl;
        }
    }

    return 0;
}