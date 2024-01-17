#ifndef WINDOWS_H
#define WINDOWS_H

#include "projetil.h"

typedef int (*compar_d_fn_t)(const void *, const void *, void *);

struct wrapper {
    compar_d_fn_t func;
    void *arg;
};

// Variável global para armazenar a estrutura wrapper durante a chamada de qsort
static struct wrapper *global_wrapper;

// Modificada para se adequar à assinatura esperada por qsort
static int compare_wrapper(const void *a, const void *b) {
    struct wrapper *wrap = global_wrapper; // Utilizando a variável global
    return wrap->func(a, b, wrap->arg);
}

void qsort_r(void *base, size_t nmemb, size_t size, compar_d_fn_t compar, void *arg) {
    struct wrapper wrap;
    wrap.func = compar;
    wrap.arg = arg;
    
    global_wrapper = &wrap; // Define a variável global antes de chamar qsort

    qsort(base, nmemb, size, compare_wrapper); // Chamada para qsort padrão
}

#endif // WINDOWS_H