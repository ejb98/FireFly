#include <stdio.h>
#include <stdlib.h>

int check_null(const char *caller_name, const char *ptr_name, void *ptr) {
    if (ptr == NULL) {
        fprintf(stderr, "%s: failed to allocate memory for field %s", caller_name, ptr_name);

        return 1;
    }

    return 0;
}