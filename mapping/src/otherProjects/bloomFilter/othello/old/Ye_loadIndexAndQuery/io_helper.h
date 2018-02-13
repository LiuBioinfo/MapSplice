// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
#ifndef __IO_HELPER_H
#define __IO_HELPER_H

template <typename keyT, typename valueT>
bool lineToKVpair(const char *s, keyT * k, int kmerLength)
{
    int tmpLength = 0;
    switch (*s) 
    {
        tmpLength ++;
        case 'A':
        case 'T':
        case 'G':
        case 'C':
            keyT ret = 0;
            while (*s == 'A' || *s == 'C' || *s =='T' || *s =='G')
            {
                ret <<=2;
                switch (*s) 
                {
                    case 'G':
                        ret ++;
                    case 'T':
                        ret ++;
                    case 'C':
                        ret ++;
                }
                s++;
                tmpLength ++;
                if(tmpLength == kmerLength)
                {
                    *k = ret;
                    return true;
                }
            }
            *k = ret;
            return true;
    }
    return false;
}
#endif