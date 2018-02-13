// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef __MUL_OTH_H
#define __MUL_OTH_H
#include "othelloV7.h"
#include "io_helper.h"
#include <cstdio>
#include <cstdlib>

using namespace std;
#define split(k,plow,phigh,spl) ((plow)=(k & ((1ULL << spl)-1)), (phigh)=(k >> spl))
template <typename keyType, typename valueType> 
class MulOth {
    vector<Othello<keyType,valueType> > vOths;
    unsigned char split;
    bool buildsucc;
public: 
    MulOth(const char * fname, unsigned char _split) {
        printf("Building Multi Othello from file %s \n",fname);
        FILE *pFile;
        pFile = fopen (fname,"r");
        vOths.clear();
        split = _split;
        for (uint32_t pl = 0; pl < (1<<_split); pl++) {
            rewind(pFile);
            printf("Reading file for keys ending with %02x/%02x\n", pl,(1<<split)-1);
            vector<keyType> keys;
            vector<valueType> values;
            while (true) {
                char buf[1024];
                if (fgets(buf,1024,pFile)==NULL) break;
                keyType k;
                valueType v;
                if (!lineToKVpair(buf, &k, &v)) break;
                uint32_t plt;
                uint64_t ph;
                split(k,plt,ph,split);
                if (plt != pl) continue;
                keys.push_back(ph);
                values.push_back(v);
            }
            printf("%d %d\n", keys.size(), values.size()); 
            Othello<keyType,valueType> oth(&keys[0], &values[0], keys.size());
            if (!oth.builtsucc) {
                printf("Build Halt!\n");
                return;
            }
            vOths.push_back(oth);
        }
        buildsucc = true;

    }
    MulOth(keyType * _keys, valueType * _values, uint64_t keycount, unsigned char _split) {
        buildsucc = false;
        vOths.clear();
        uint32_t plmax = (1<<_split) - 1;
        split = _split;
        for (uint32_t pl = 0; pl <= plmax; pl++) {
            vector<keyType> keys;
            vector<valueType> values;
            for (uint64_t i = 0 ; i < keycount; i++) {
                uint32_t plt;
                uint64_t ph;
                split(_keys[i],plt,ph,split);
                if (plt != pl) continue;
                keys.push_back(ph);
                values.push_back(_values[i]);
            }
            Othello<keyType,valueType> oth(&keys[0], &values[0], keys.size());
            if (!oth.builtsucc) {
                printf("Build Halt!\n");
                return;
            }
            vOths.push_back(oth);
            //oth.printall();
        }
        buildsucc = true;
    }
    inline valueType query(const keyType &k) {
        uint32_t pl;
        uint64_t ph;
        split(k,pl,ph,split);
        return vOths[pl].query(ph);
    }
    void printall () {
        printf("Printall ...\n");
        for (auto V:vOths)
            V.printall();
    }
    void writeToFile(const char* fname) {
        FILE *pFile;
        pFile = fopen (fname, "wb");
        unsigned char buf0x20[0x20];
        memset(buf0x20,0, sizeof(buf0x20));
        uint32_t versionID = VERSIONID;
        memcpy(buf0x20,&versionID, sizeof(uint32_t));
        uint32_t split32 = split;
        memcpy(buf0x20+0x4, & split32, sizeof(uint32_t));
        fwrite(buf0x20,sizeof(buf0x20),1,pFile);
        for (int i = 0 ; i <(1<<split); i++) {
            vOths[i].exportInfo(buf0x20);
            fwrite(buf0x20,sizeof(buf0x20),1,pFile);
        }
        for (int i = 0 ; i <(1<<split); i++) {
            fwrite(vOths[i].m,sizeof(valueType),vOths[i].othelloHashUpperBound,pFile);
        }
        fclose(pFile);
    }
     
    MulOth(const char* fname) {
        buildsucc = false;
        printf("Read from binary file %s\n", fname);
        FILE *pFile;
        pFile = fopen (fname, "rb");
        uint32_t versionID = VERSIONID;
        uint32_t compversion;
        unsigned char buf0x20[0x20];
        fread(buf0x20,sizeof(buf0x20),1,pFile);
#ifndef NO_VERSION_CHECK        
        if (memcmp(buf0x20,&versionID,sizeof(uint32_t))) {
            printf("Wrong version number\n");
            fclose(pFile);
            return;
        }
#endif        
        uint32_t split32;
        memcpy(&split32, buf0x20+0x4 , sizeof(uint32_t));
        split = split32;
        for (int i = 0 ; i < (1<<split); i++) {
           fread(buf0x20,sizeof(buf0x20),1,pFile);
           Othello<keyType,valueType> *oth;
           oth = new Othello<keyType,valueType>(buf0x20);
           vOths.push_back(*oth);
        }
        for (int i = 0 ; i < (1<< split); i++) {
           fread(vOths[i].m,sizeof(valueType), vOths[i].othelloHashUpperBound, pFile);
        }
        fclose(pFile);
        buildsucc = true;
    } 
};


//MulOthello binary file descriptor
//0x00 : uint32_t signature MulOth version 
//0x04 : uint32_t splitbit
//0x20 : OthInfo1
//0x40 : OthInfo2
//...
//offset1 : Oth[0].m
//offset2 = offset1 + Oth[0].hashupperbound = Oth[1].m
//...

//OthInfo: 32 Byte
//0x00 : uint64_t hash1
//0x08 : uint64_t hash2
//0x10 : uint32_t mask1
//0x14 : uint32_t mask2
//0x18 : uint64_t m.offset

#endif
