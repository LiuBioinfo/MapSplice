/*
 * bloom.h
 *
 *  Created on: 2012-2-22
 *      Author: xiaojay
 */

#ifndef BLOOM_H_
#define BLOOM_H_
#include <vector>
#include "hashFunc.h"

class Bloom
{
public:
    Bloom(int size , std::vector<HashFun*> hashfunclist)
    {
        assert(hashfunclist.size()>0);
        this->size = size;
        this->hashfunclist = hashfunclist;
        this->arr = new char [size];
    }

    ~Bloom()
    {
        if(this->arr!=NULL)
        {
            delete this->arr;
        }
    }

    void add(const char * text)
    {
        int nfunc = hashfunclist.size();
        long code = 0;
        for(int i=0;i<nfunc;i++)
        {
            code = hashfunclist.at(i)->gethashval(text);

            if(code/CHARBITSIZE>size) return;
            else
            {
                setbit(code);
            }
        }
    }

    bool check(const char * text)
    {
        int nfunc = hashfunclist.size();
        long code = 0;
        for (int i=0;i<nfunc;i++)
        {
            code = hashfunclist.at(i)->gethashval(text);
            if(code/CHARBITSIZE>size) return false;
            else
            {
                if(getbit(code)) 
                    continue;
                else
                    return false;
            }
        }
        return true;
    }

private:
    const static int CHARBITSIZE = 8;
    int size;
    char * arr;
    std::vector<HashFun*> hashfunclist;

    inline void setbit(long code)
    {
        arr[code/CHARBITSIZE] |= (1<<(code%CHARBITSIZE));
    }

    inline bool getbit(long code)
    {
        if(!(arr[code/CHARBITSIZE] & (1<<(code%CHARBITSIZE))))
        {
            return false;
        }
        return true;
    }
};
#endif