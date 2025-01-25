

#ifndef __SEMINMAX_H__
#error Include only from SEminmax.h
#endif

#ifndef minmax_type
#error minmax_type should be defined
#endif
#ifndef minmax_name
#error minmax_name should be defined
#endif


SE_inline minmax_type minmax_name(min)(minmax_type a,
                                                    minmax_type b)
{
    if(a <= b ) return a;
    else        return b;
}

SE_inline minmax_type minmax_name(max)(minmax_type a,
                                                    minmax_type b)
{
    if(a >= b ) return a;
    else        return b;
}

SE_inline minmax_type minmax_name(min3)(minmax_type a,
                                                     minmax_type b,
                                                     minmax_type c)
{
    minmax_type r;

    if(a < b) r = a;
    else      r = b;

    if(c < r) r = c;

    return r;
}

SE_inline minmax_type minmax_name(max3)(minmax_type a,
                                                     minmax_type b,
                                                     minmax_type c)
{
    minmax_type r;

    if(a > b) r = a;
    else      r = b;

    if(c > r) r = c;

    return r;
}

SE_inline minmax_type minmax_name(min4)(minmax_type a,
                                                     minmax_type b,
                                                     minmax_type c,
                                                     minmax_type d)
{
    minmax_type r;

    if(a < b) r = a;
    else      r = b;

    if(c < r) r = c;

    if(d < r) r = d;

    return r;
}

SE_inline minmax_type minmax_name(max4)(minmax_type a,
                                                     minmax_type b,
                                                     minmax_type c,
                                                     minmax_type d)
{
    minmax_type r;

    if(a > b) r = a;
    else      r = b;

    if(c > r) r = c;

    if(d > r) r = d;

    return r;
}
