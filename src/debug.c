#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>

void print_trace(){
    void *array[20];
    size_t size;
    char **strings;
    int i;
    size = backtrace(array,20);
    strings = backtrace_symbols(array,size);
    printf("Obtained %zd stack frames.\n",size);
    for(i=0;i<size;i++)
        printf("%s\n",strings[i]);
    free(strings);
}
