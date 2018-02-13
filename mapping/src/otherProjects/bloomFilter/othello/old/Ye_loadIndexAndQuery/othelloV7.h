#ifndef __OTHELLO_V3_H
#define __OTHELLO_V3_H
//Supports Othello remote update
#include <cstdint>
#include <ctime>
#include <cstring>
#include <algorithm>
#include <vector>
#include <map>
//#include "hashes.h"
#define COUNTERLENGTH 6
#define COUNTER_ADD_BITS 3
#define COUNTER_ADD_MESH 0x7
#define COUNTER_ALL_MESH 0x1FF
#include <vector>
#define MAX_REHASH 10
using namespace std;
template<typename keyType>
void inline CRC32CASM1(uint32_t s1, keyType * k0, uint32_t &ret1, int hashshr, int mask) {
//    static int mask = ( 1 << (hashlength) ) -1;
    static_assert(
        sizeof(keyType)==4
        || sizeof(keyType)==8
        || sizeof(keyType)==12
        || sizeof(keyType)==16
        || sizeof(keyType)==20
        || sizeof(keyType)==24
        || sizeof(keyType)==28
        || sizeof(keyType)==32
        , "keyType length should be 32/64/96/128/160/192/224/256 bits ");
    uint32_t crc1 = ~0;
    uint64_t *k;
    k = (uint64_t* ) k0;
    if (sizeof(keyType)>=8)
        asm(".byte 0xf2, 0x48, 0xf, 0x38, 0xf1, 0xf1;"  :"=S"(crc1)   :"0"(crc1), "c" ((*k)+s1)  );

    if (sizeof(keyType)>=16) {
        s1 = ((((uint64_t) s1) * s1 >> 16) ^ (s1 << 2));
        k++;
        asm(".byte 0xf2, 0x48, 0xf, 0x38, 0xf1, 0xf1;"  :"=S"(crc1)   :"0"(crc1), "c" ((*k)+s1)  );
    }
    if (sizeof(keyType)>=24) {
        s1 = ((((uint64_t) s1) * s1 >> 16) ^ (s1 << 2));
        k++;
        asm(".byte 0xf2, 0x48, 0xf, 0x38, 0xf1, 0xf1;"  :"=S"(crc1)   :"0"(crc1), "c" ((*k)+s1)  );
    }
    if (sizeof(keyType)>=32) {
        s1 = ((((uint64_t) s1) * s1 >> 16) ^ (s1 << 2));
        k++;
        asm(".byte 0xf2, 0x48, 0xf, 0x38, 0xf1, 0xf1;"  :"=S"(crc1)   :"0"(crc1), "c" ((*k)+s1)  );
    }
    if (sizeof(keyType) & 4) {
        uint64_t k32;
        k32 =  ((*k)) & 0xFFFFFFFFULL;
        asm( ".byte 0xf2, 0xf, 0x38, 0xf1, 0xf1;"  :"=S"(crc1)    :"0" (crc1),"c" ((k32)+s1)  );

    }

    asm( ".byte 0xf2, 0xf, 0x38, 0xf1, 0xf1;"  :"=S"(crc1)    :"0" (crc1),"c" (s1)  );
//    crc1 ^= (crc1 >> (HASHLENGTH ^ (7&s1)));
    crc1 ^= (crc1 >> (hashshr));
    if (sizeof(keyType)==4)
        ret1 = mask & (crc1 ^ ((uint32_t) *k));
    else
        ret1 = mask & (crc1 ^ (*k >> 32) ^ ((uint32_t) *k));
    return;
}



template<typename keyType, typename valueType>
class Othello {
public:
    valueType * m;
#ifdef ENABLE_UPDATE    
    uint8_t * Z1;
    uint8_t * Z2;
#endif    
    uint32_t s1, s2;
    uint8_t *filled;
    uint32_t hashshr1, mask1;
    uint32_t hashshr2, mask2, m2shift;
private:
    keyType *keys;
    valueType *values;
    uint32_t *next1;    //next1[i]: h1(keys[i]) = h1(keys[next[i]]);
    uint32_t *next2;    //next2[i]: h1(keys[i]) = h1(keys[next[i]]);
    uint32_t *first;    //h1(keys[first[i]]) =h1(keys[last[i]]) = i // if i < half     //h2(keys[first[i]]) = i // if i >=half
    uint32_t *fa;

public:
    uint32_t othelloHashUpperBound;
    uint32_t keycount;
    uint32_t keycountPreserve;
#ifndef ENABLE_UPDATE
    void exportInfo(unsigned char * v) {
         memset(v,0,0x20);
         memcpy(v,&s1,sizeof(uint32_t));
         memcpy(v+8,&s2,sizeof(uint32_t));
         memcpy(v+0x10,&mask1, sizeof(uint32_t));
         memcpy(v+0x14,&mask2,sizeof(uint32_t));
    }
    Othello(unsigned char *v) {
         memcpy(&s1,v,sizeof(uint32_t));
         memcpy(&s2,v+8,sizeof(uint32_t));
         memcpy(&mask1, v+0x10, sizeof(uint32_t));
         memcpy(&mask2, v+0x14, sizeof(uint32_t));
         hashshr1 = s1 & 7;
         hashshr2 = s2 & 7;
         m2shift = mask1+1;
         othelloHashUpperBound = mask1+mask2+2;
         m= (valueType*) valloc(othelloHashUpperBound*sizeof(valueType));
    }
#endif    
private:
    void findnewhash() {
        s1 = rand();
        s2 = rand();
        hashshr1 = s1 & 7;
        hashshr2 = s2 & 7;
        memset(m,0,sizeof(valueType)* othelloHashUpperBound);
#ifdef ENABLE_UPDATE        
        memset(Z1,0,1 << COUNTERLENGTH);
        memset(Z2,0,1 << COUNTERLENGTH);
#endif        
        memset(filled, 0, keycountPreserve>>3);
        trycount++;
    }
    uint32_t getfa(uint32_t k) {
        if (fa[k] & 0xF0000000) return fa[k] = k;
        if (fa[k]!=k) return fa[k] = getfa(fa[k]);
        else return k;
    }
    void merge(uint32_t k1, uint32_t k2) {
        fa[getfa(k1)] = getfa(k2);
    }
#ifdef ENABLE_UPDATE    
    void growkeycountPreserve() {
        uint32_t oct = keycountPreserve;
        keycountPreserve *= 2;
        printf("Grow Preserve %d %d\n", oct, keycountPreserve);
        uint8_t * newfilled = (uint8_t *) valloc((keycountPreserve >> 3)+1);
        uint32_t * newnext1 = (uint32_t *) valloc(keycountPreserve * sizeof(uint32_t ));
        uint32_t * newnext2 = (uint32_t *) valloc(keycountPreserve * sizeof(uint32_t ));

        keyType * newkeys = (keyType * ) valloc(keycountPreserve * sizeof(keyType));
        valueType * newvalues = (valueType * ) valloc(keycountPreserve * sizeof(valueType));
        memcpy(newfilled,filled,(oct>>3)+1);
        memcpy(newnext1,next1,sizeof(uint32_t)*oct);
        memcpy(newnext2,next2,sizeof(uint32_t)*oct);
        memcpy(newkeys,keys,sizeof(keyType)*oct);
        memcpy(newvalues,values,sizeof(valueType)*oct);

        free(filled);
        filled=newfilled;
        free(next1);
        next1=newnext1;
        free(next2);
        next2=newnext2;
        free(keys);
        keys=newkeys;
        free(values);
        values = newvalues;
    }
#endif    
    bool trybuild(keyType *keys, valueType *values, uint32_t keycount);
    void dfs(uint32_t);
    struct timespec tm1,tm2;
    void filldfs(uint32_t, valueType, uint32_t);
#ifdef ENABLE_UPDATE    
    vector< pair< uint32_t, valueType> > updateBuffer;
    void clearUpdateBuffer();
    void modifyValueToDesire(keyType k, valueType v);
#ifdef ENABLE_FRIENDS    
    void notifyfriends();
#endif    
    bool keyexists(keyType k, bool & addsucc, valueType v, bool tryadd);
#endif    

public:
    Othello(keyType * _keys, valueType *_values, uint32_t keycount);
    Othello(valueType * _m, uint8_t *Z1, uint8_t *Z2, uint32_t _s1, uint32_t _s2, uint32_t _mask1, uint32_t _mask2);
    void printkeys() {
        for (int i = 0 ; i < keycount; i++) {
            uint32_t h1,h2;
            keyType k = keys[i];
            get_hash(k,h1,h2);
            printf("%08llx:%d %x-%x==>%d.%d; q%x  %x ^ %x\n", keys[i],values[i],h1,h2,getfa(h1),getfa(h2), query(keys[i]), m[h1],m[h2]);
        }
        printf("===========\n");
    }
    void printall() {
        printf("%08x %08x %08x %08x\n", mask1,mask2, s1,s2);
        for (int i = 0 ; i < 1+mask1+1+mask2; i++) {
            printf("%02x",m[i]);
            if ((i+1)%4==0) printf(" ");
            if ((i+1)%16==0) printf("\n");
        }
       printf("\n"); 
    }
    void inline get_hash_1(const keyType &v, uint32_t &ret1) {
        CRC32CASM1(s1,&v,ret1,hashshr1, mask1);
    }
    void inline get_hash_2(const keyType &v, uint32_t &ret2) {
        CRC32CASM1(s2,&v,ret2,hashshr2,mask2);
        ret2 += m2shift;
        //printf("%llx %x \n", v, ret2);
    }
    void inline get_hash(const keyType &v, uint32_t &ret1, uint32_t &ret2) {
        get_hash_1(v,ret1);
        get_hash_2(v,ret2);
        //printf("%llx: %llx %x %x : %x %x \n",((uint64_t) this), v, ret1, ret2,s1,s2);
    }

    inline valueType query(const keyType & k);
#ifdef ENABLE_UPDATE    
private:
    bool update(keyType &k, valueType &v);
public:
    bool batchupdate(keyType *k, valueType *v, uint32_t cnt);
    bool batchupdateOrInsert(keyType *k, valueType *v, uint32_t cnt, bool skipNewkeyCheck);
#endif
    bool builtsucc;
    uint32_t trycount;
#ifdef ENABLE_SUSCEPTBILITY    
    double getSusceptibility();
#endif
#ifdef ENABLE_FRIENDS    
    vector< Othello<keyType,valueType> * > friends;
    vector< void *> friendm;
    vector< pair< void *, void *> > friendZ;
    void addfriend( Othello<keyType,valueType> * f);
#endif    
};

#ifdef ENABLE_SUSCEPTBILITY
template<typename keyType, typename valueType>
double Othello<keyType,valueType>::getSusceptibility() {
    map<int,int> VK;
    for (int i = 0 ; i < othelloHashUpperBound; i++)
        VK[getfa(i)] ++;
    double ans = 0;
    for (auto v : VK) {
        ans += (v.second * v.second);
    }
    return ans/othelloHashUpperBound;
}

#endif 

template<typename keyType, typename valueType>
void Othello<keyType,valueType>::filldfs(uint32_t root, valueType flip, uint32_t grandpa) {
    uint32_t v2 = first[root];
    m[root] ^= flip;
//    printf("dfs %d ,grandpa %d:\n", root,grandpa);
#ifdef ENABLE_UPDATE    
#define ADDLOG(x,y) updateBuffer.push_back(make_pair(x,y));
#else
#define ADDLOG(x,y) ;
#endif    
    ADDLOG(root,m[root]);
    while ((v2 & 0xF0000000)==0) {
        uint32_t h1, hx;
        keyType kk = keys[v2];
        get_hash(kk,h1,hx);
        if (h1 != grandpa) {
            //      	printf("grandpa%d root%d key %d %lx : %d - %d : \n ", grandpa, root, v2,kk,h1,hx );
            m[h1] ^= flip;
            ADDLOG(h1,m[h1]);
            uint32_t v1 = first[h1];
            while ((v1 & 0xF0000000) ==0 ) {
                if (v1 != v2) {
                    keyType kk = keys[v1];
                    uint32_t h2,hx;
                    get_hash(kk,hx,h2);
                    //                printf("grandpa %d,root%d,h1 %d  ==>key %d %lx : %d - %d : \n", grandpa,root, h1, v1,kk,hx,h2 );
                    filldfs(h2, flip, h1);
                }
                v1 = next1[v1];
            }
        }
        v2 = next2[v2];
    }
}

#ifdef ENABLE_UPDATE
template<typename keyType, typename valueType>
bool Othello<keyType, valueType>::keyexists(keyType k, bool &addsucc , valueType v, bool tryadd) {
    //exzamine whether k exists in the graph, if not,
    //try to add a key into graph, but not modify the values,
    //this key must be unique new key.
    //return false if the key does not belong to keyset before such test.
    //return false and, but addsucc=false if will introduce a circle.
    //TODO

    uint32_t h1,h2;
    get_hash(k,h1,h2);
    uint32_t f1,f2;
    f1 = getfa(h1);
    f2 = getfa(h2);
    if (f1!=f2) {
        // this is a new key!
        if (addsucc = tryadd ) {
            if (keycount *2 >= keycountPreserve) {
                growkeycountPreserve();
            }
//			printf("%d %d\n",keycountPreserve,keycount);

            next1[keycount] = first[h1];
            next2[keycount] = first[h2];
            first[h2] = first[h1] = keycount;
            merge(h1,h2);
            keys[keycount] = k;
            values[keycount] = v;
            keycount ++;

        }

        return false;
    }
    //f1 == f2 : this is an old key, or  a new key but hash value (h1,h2) exists.
    uint32_t id = first[h1];
    while ((id & 0xF0000000) == 0) {
        if (keys[id] == k) {
            addsucc=false;
            return true;
        }
        uint32_t hv2;
        get_hash_2(keys[id], hv2);
        if (hv2 == h2) {
            addsucc = false;
            return false;
        }
        id = next1[id];
    }
    //if goes here, means we will have a circle;
    return false;
}
template<typename keyType, typename valueType>
void Othello<keyType,valueType>::clearUpdateBuffer() {
    updateBuffer.clear();
}

template<typename keyType, typename valueType>
void Othello<keyType,valueType>::modifyValueToDesire( keyType k, valueType v) {
    //TODO
    uint32_t h1,h2;
    get_hash(k,h1,h2);
#ifdef TEST_MODIFY
//    printf("Modify %llx: %x %x \n",k,m[h1]^m[h2],v);
#endif
    filldfs(h2, m[h2]^m[h1]^v, h1);
}
#endif 

#ifdef ENABLE_FRIENDS
template<typename keyType, typename valueType>
void Othello<keyType,valueType>::notifyfriends() {
    //this can only notify without S1S2 change.
    for (int i = 0 ; i< friendm.size(); i++) {
//		printf("Notify %llx of %d changes", (uint64_t ) friendm[i], updateBuffer.size());
        valueType * mother = (valueType *) friendm[i];
        uint8_t * Z1other = (valueType *) friendZ[i].first;
        uint8_t * Z2other = (valueType *) friendZ[i].second;
        for (auto p: updateBuffer) {
            uint32_t id = p.first & COUNTER_ALL_MESH;
            Z1other[id>>COUNTER_ADD_BITS] ^= (1<<(id & COUNTER_ADD_MESH));
        }
        __asm__ __volatile__ ("" ::: "memory");
        for (auto p: updateBuffer) {
            mother[p.first] = p.second;
        }
        __asm__ __volatile__ ("" ::: "memory");
        for (auto p: updateBuffer) {
            uint32_t id = p.first & COUNTER_ALL_MESH;
            Z2other[id>>COUNTER_ADD_BITS] ^= (1<<(id & COUNTER_ADD_MESH));
        }
    }

}

template<typename keyType, typename valueType>
void Othello<keyType,valueType>::addfriend( Othello<keyType,valueType> *f) {
    void * fm = (void *) f->m;
    if (find(friendm.begin(),friendm.end(),fm)==friendm.end()) {
        friendm.push_back(fm);
        friendZ.push_back(make_pair((void *) f->Z1, (void *) f->Z2));
    }
    friends.push_back(f);
}

#endif
#ifdef ENABLE_UPDATE
template<typename keyType, typename valueType>
bool Othello<keyType, valueType>::batchupdate(keyType *k, valueType *v, uint32_t cnt) {
    //return true if all these cnt keys are successfully updated. if not, this is to say adding a new key will introduce a cycle in the graph. return false. In this case, should consturct a new Othello instance.
    //Note: this function should only be called when it is sure that the keys already exists in keyset.
    clearUpdateBuffer();
    for (int i = 0 ; i < cnt; i++) {
        bool addsucc;
#ifdef TEST_MODIFY
//            printf("::%d:: "            ,i);
#endif
        bool keyexist = keyexists(k[i], addsucc,v[i],true);
        if ((!keyexist) && !(addsucc)) {
            //TODO: clear update history.
            return false;
        }
        modifyValueToDesire(k[i],v[i]);

    }
    notifyfriends();
    return true;
}

template<typename keyType, typename valueType>
bool Othello<keyType, valueType>::batchupdateOrInsert(keyType *k, valueType *v, uint32_t cnt, bool skipNewkeyCheck) {
    // returns true if this updateorinsert does not include rehash.
    int bkkeycount = keycount;
    if (batchupdate(k,v,cnt)) {

        //add keys to keyset log.
        return true;
    }
    keycount = bkkeycount;
    bool tmpbool;
    valueType tmpv;
    for (int i = 0; i < cnt; i++) {
        bool addthiskey = true;
        if (!skipNewkeyCheck)
            addthiskey = !keyexists(k[i],tmpbool,tmpv,false);
        if (addthiskey) {
            if (keycount * 2>= keycountPreserve) {
                growkeycountPreserve();
            }
            keys[keycount] = k[i];
            values[keycount] = v[i];
            keycount ++;
        }
    }

    findnewhash();
    builtsucc = false;
    trycount = 0;

#ifdef PREFER_KILO
    printf("Re=Building Othello %dk Keys, %dk/%dk nodes , keysize %d, valuesize %d \n",keycount >> 10, (mask1+1)>>10, (mask2+1)>>10, sizeof(keyType), sizeof(valueType));
#else
    printf("Re=Building Othello %d Keys, %d/%d nodes, keysize %d, valuesize %d \n",keycount, mask1+1,mask2+1, sizeof(keyType),sizeof(valueType) );
#endif
    while ((!builtsucc)&&(trycount < MAX_REHASH)) {
        builtsucc = trybuild(keys,values, keycount);
        if (builtsucc) break;
        trycount ++;
        findnewhash();
    }
    printf("S1 %x  S2 %x, succ %d after %d tries\n", s1,s2,builtsucc, trycount);



    if (!builtsucc) {
        printf("!!!! not built succ!");
    }
    return false;
}
#endif 




template<typename keyType, typename valueType>
valueType Othello<keyType,valueType>::query(const keyType &k) {
    uint32_t h1, h2;
#ifndef PREFETCH
    get_hash(k,h1,h2);
#else
    get_hash_1(k,h1);
#endif
    uint8_t v1,v2;
    uint8_t c1 = h1 & COUNTER_ALL_MESH;
    while (true) {
#ifdef ENABLE_UPDATE
        uint8_t z1 = Z1[c1>>COUNTER_ADD_BITS]>>(c1 & COUNTER_ADD_MESH);
#endif
        __asm__ __volatile__ ("" ::: "memory");
        v1 = (m[h1]); // >> ((h1 & 0x7)<<VALUE_ADD_BITS));
#ifdef PREFETCH
        get_hash_2(k,h2);
#endif
        v2 = (m[h2]);// >> ((h2 & 0x7)<<VALUE_ADD_BITS));
        __asm__ __volatile__ ("" ::: "memory");
//    uint8_t c2 = h2 & 0xFFF;
#ifdef ENABLE_UPDATE       
        uint8_t z2 = Z2[c1>>COUNTER_ADD_BITS]>>(c1 & COUNTER_ADD_MESH);
#endif
#ifdef PRINT_BEFORE_QUERY
        printf("%x %x--> %x %x \n",h1,h2,m[h1],m[h2]);
#endif      
#ifdef ENABLE_UPDATE        
        if (((z1^z2)&1)==0) 
#endif            
            return v1^v2;
    };
}


template<typename keyType, typename valueType>
Othello<keyType,valueType>::Othello(keyType * _keys, valueType *_values, uint32_t _keycount) {
    int hl1 = 0;
    int hl2 = 0 ;
    while ((1<<hl2) <  _keycount * 1) hl2++;
    while ((1<<hl1) < _keycount* 1.333334) hl1++;
    mask1 = (1<<hl1) -1;
    mask2 = (1<<hl2) -1;
    hashshr1 = hl1 ^ ( 7 & s1);
    hashshr2 = hl2 ^ ( 7 & s2);
    m2shift = mask1+1;

    keycount = _keycount;
//    keycountPreserve = (1<<HASHLENGTH) * 3 / 2;
#ifdef ENABLE_UPDATE
    keycountPreserve = keycount * 2;
#else
    keycountPreserve = keycount;
#endif    
    othelloHashUpperBound = m2shift + mask2 + 1;
    tm1.tv_sec = 0;
    tm1.tv_nsec = 5;
    m = (valueType *) valloc(othelloHashUpperBound* sizeof(valueType));
#ifdef ENABLE_UPDATE
    Z1 = (uint8_t *) valloc(sizeof(uint8_t)  << COUNTERLENGTH);
    Z2 = (uint8_t *) valloc(sizeof(uint8_t)  << COUNTERLENGTH);
#endif    
    filled = (uint8_t *) valloc((keycountPreserve >> 3)+1);
    next1 = (uint32_t *) valloc(keycountPreserve * sizeof(uint32_t ));
    next2 = (uint32_t *) valloc(keycountPreserve * sizeof(uint32_t ));

    first = (uint32_t *) valloc(othelloHashUpperBound * sizeof(uint32_t ));
    keys = (keyType * ) valloc(keycountPreserve * sizeof(keyType));
    values = (valueType * ) valloc(keycountPreserve * sizeof(valueType));
    fa = (uint32_t * ) valloc(othelloHashUpperBound*sizeof(uint32_t));
    memcpy(keys, _keys, sizeof(keyType)*keycount);
    memcpy(values, _values, sizeof(valueType)*keycount);
    builtsucc = false;
    trycount = 0;
    findnewhash();
#ifdef PREFER_KILO
    printf("Building Othello %dk Keys, %dk / %dk nodes, keysize %d, valuesize %d \t",keycount >> 10, (mask1+1)>>10, (mask2+1)>>10, sizeof(keyType), sizeof(valueType));
#else
    printf("Building %d Keys, %d/%d nodes, Ksize %dB, Vsize %dB\t",keycount, mask1+1, mask2+1, sizeof(keyType),sizeof(valueType) );
#endif
    while ((!builtsucc)&&(trycount < MAX_REHASH)) {
        builtsucc = trybuild(keys,_values, keycount);
        if (builtsucc) break;
        trycount ++;
        findnewhash();
    }
#ifdef PRINT_BEFORE_QUERY
    printkeys();
#endif
    free(first);
    free(fa);
    free(next1);
    free(next2);
    free(keys);
    free(values);
    printf("S1 %x S2 %x, succ %d after %d tries\n", s1,s2,builtsucc, trycount);
}

#ifdef ENABLE_FRIENDS
template<typename keyType, typename valueType>
Othello<keyType,valueType>::Othello(valueType * _m, uint8_t * _Z1, uint8_t * _Z2, uint32_t _s1, uint32_t _s2, uint32_t _mask1, uint32_t _mask2) {
    s1 = _s1;
    s2 = _s2;
    mask1 = _mask1;
    mask2 = _mask2;
    m2shift = mask1+1;
    int hl1 = 0;
    while (_mask1) {
        hl1++;
        _mask1 >>=1;
    };
    int hl2 = 0;
    while (_mask1) {
        hl2++;
        _mask2 >>=1 ;
    };
    hashshr1 = hl1 ^ ( 7 & s1);
    hashshr2 = hl2 ^ ( 7 & s2);
    builtsucc = true;
    m = _m;
    Z1 =  _Z1;
    Z2 = _Z2;
    tm1.tv_sec = 0;
    tm1.tv_nsec = 0;
    othelloHashUpperBound = m2shift + mask2 + 1;
    keys = NULL;
}

#endif 


template<typename keyType, typename valueType>
bool Othello<keyType,valueType>::trybuild(keyType *keys, valueType *values, uint32_t keycount) {
    memset(fa, 0xFF, sizeof(uint32_t)*othelloHashUpperBound);
    memset(next1,0xFF,sizeof(uint32_t)* keycount);
    memset(next2,0xFF,sizeof(uint32_t)* keycount);
    memset(first,0xFF,sizeof(uint32_t)* othelloHashUpperBound);
    uint32_t h1, h2;
    for (uint32_t i = 0; i < keycount; i++) {
        keyType kk = keys[i] ;
        get_hash(kk, h1, h2);
#ifdef TEST
        printf("%d %llx %x %x :%x %x\n",i,kk,h1,h2,s1,s2);
#endif
        if (getfa(h1) == getfa(h2)) return false;
        merge(h1,h2);
        if (first[h1] & 0xF0000000)
            first[h1] = i;
        else
        {
            next1[i] = first[h1];
            first[h1] = i;
        }
        if (first[h2] & 0xF0000000)
            first[h2] = i;
        else
        {
            next2[i] = first[h2];
            first[h2] = i;
        }
    }
    memset(filled,0,(keycountPreserve >> 3) +1);
    for (uint32_t i = 0; i < keycount; i++) {
        bool isfilled = filled[i>>3] & (1U << (i & 0x7));
        if (isfilled) continue;
        uint32_t h1,h2;
        keyType kk = keys[i] ;
        get_hash(kk, h1,h2);
        dfs(h1);
    }
#ifdef TEST
    for (uint32_t i = 0 ; i < keycount; i++) {
        uint32_t h1,h2;
        keyType kk = keys[i] ;
        get_hash(kk, h1,h2);
        printf("%4d %16llx %8x %8x %4x %4x %4x\n",i, keys[i], h1,h2,
               (uint8_t)  (m[h1]), (uint8_t) (m[h2]),
               (uint8_t) ((m[h1]) ^ (m[h2])));
    }
#endif
    return true;
}

template<typename keyType, typename valueType>
void Othello<keyType,valueType>::dfs(uint32_t h) {
    uint32_t i = first[h];
    while ( (i & 0xF0000000) == 0) {
        keyType p;
        p = keys[i] ;// & KEYMESH;
        uint32_t hx,h2;
        get_hash(p,hx,h2);
        valueType v;
        v = values[i];
        uint8_t filledbit = 1U << (i & 0x7);
        if (filled[i>>3] & filledbit) {
            i = next1[i];
            continue;
        }
        v ^= m[h];
        m[h2] = v;
        filled[i>>3] |= filledbit;
        uint32_t j = first[h2];
        while ((j & 0xF0000000)==0) {
            if (j==i) {
                j = next2[j];
                continue;
            }
            keyType q;
            q = keys[j];
            uint32_t h1,hy;
            get_hash(q,h1,hy);
            valueType u;
            u = values[j];
            u ^= m[h2];
            m[h1] = u;
            filled[j>>3] |= (1U << (j&0x7));
            dfs(h1);
            j = next2[j];
        }
        i = next1[i];
    }
}


#endif

