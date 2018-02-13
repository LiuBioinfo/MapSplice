#ifndef __IO_HELPER_H
#define __IO_HELPER_H

template <typename keyT, typename valueT>
bool lineToKVpair(char *s, keyT * k, valueT *v) {
    switch (*s) {
    case 'A':
    case 'T':
    case 'G':
    case  'C' :
        keyT ret = 0;
        while (*s == 'A' || *s == 'C' || *s =='T' || *s =='G') {
            ret <<=2;
            switch (*s) {
            case 'G':
                ret ++;
            case 'T':
                ret ++;
            case 'C':
                ret ++;
            }
            s++;
        }
        *k = ret;
        valueT tv;
        sscanf(s,"%d",&tv);
        *v = tv;
        return true;

    }
    return false;
}

#endif

