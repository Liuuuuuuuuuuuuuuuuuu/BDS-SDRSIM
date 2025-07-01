#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

int main(int argc, char *argv[]) {
    const char *filename = "beidou_b1i.bin";
    FILE *f = fopen(filename, "rb");
    if (!f) {
        perror("Unable to open file");
        return 1;
    }

    const size_t samples = 10000; // number of IQ pairs
    int16_t *buffer = malloc(samples * 2 * sizeof(int16_t));
    if (!buffer) {
        perror("Unable to allocate buffer");
        fclose(f);
        return 1;
    }

    size_t read = fread(buffer, sizeof(int16_t), samples * 2, f);
    fclose(f);
    if (read < samples * 2) {
        fprintf(stderr, "File contains only %zu I/Q values\n", read);
    }

    long long sum_i = 0, sum_q = 0;
    int max_abs_i = 0, max_abs_q = 0;

    for (size_t n = 0; n < read / 2; n++) {
        int16_t i = buffer[2*n];
        int16_t q = buffer[2*n + 1];
        printf("Sample %zu: I=%d, Q=%d\n", n+1, i, q);
        sum_i += i;
        sum_q += q;
        if (abs(i) > max_abs_i) max_abs_i = abs(i);
        if (abs(q) > max_abs_q) max_abs_q = abs(q);
    }

    free(buffer);

    double avg_i = (double)sum_i / (read / 2);
    double avg_q = (double)sum_q / (read / 2);

    printf("I average: %.2f, I max abs: %d\n", avg_i, max_abs_i);
    printf("Q average: %.2f, Q max abs: %d\n", avg_q, max_abs_q);

    if (fabs(avg_q) < fabs(avg_i) * 0.1 && max_abs_q < max_abs_i / 10) {
        printf("Q values are much smaller than I. Possible issue with Q channel.\n");
    } else {
        printf("I and Q have comparable magnitude. Data looks fine.\n");
    }

    return 0;
}
