#include <util.h>

int bin_coeff(int n, int r) {
    int b;

    if ((r < 0) || (n < r)) {
        return (0);
    }

    if ((2 * r) > n) {
        r = n - r;
    }

    b = 1;

    if (r > 0) {
        for (int i = 0; i <= (r - 1); i = i + 1) {
            b = (b * (n - i)) / (i + 1);
        }
    }

    return (b);
}

void print_line() {
    printf(
        "---------------------------------------------------------------------"
        "-------------\n");
}

int bisearch(int* sorted, const void* target, int size,
             int (*compare)(const void* key1, const void* key2)) {
    int left, right;
    left = 0;
    right = size - 1;

    while (left <= right) {
        int middle = (left + right) / 2;

        switch (compare((sorted + (middle)), target)) {
            case -1:
                left = middle + 1;
                break;

            case 1:
                right = middle - 1;
                break;

            case 0:
                return middle;
        }
    }

    return -1;
}

void dump_uname() {
    struct utsname uts;
    uname(&uts);
    printf("sysname: %s\n", uts.sysname);
    printf("nodename: %s\n", uts.nodename);
    printf("release: %s\n", uts.release);
    printf("version: %s\n", uts.version);
    printf("machine: %s\n", uts.machine);
}

int program_header(int ac, char** av) {
    int    val = 0;
    time_t starttime;
    int    i;
    printf(
        "#####################################################################"
        "###################################\n");
    printf("Running: ");

    for (i = 0; i < ac; ++i) {
        printf(" %s", av[i]);
    }

    printf("\n");
    dump_uname();
    (void)time(&starttime);
    printf(
        "#####################################################################"
        "###################################\n");
    return val;
}
