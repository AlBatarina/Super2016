/* Written by Dmitry Chestnykh. Public domain. */
#include <err.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <openssl/blowfish.h>
#include <omp.h>

/* Encryption and decryption */

void
cipher(int what, uint8_t *dst, uint8_t *src, uint8_t key)
{
	BF_KEY bk;

	BF_set_key(&bk, 1, &key);
	BF_ecb_encrypt(src, dst, &bk, what);
}

void
double_encrypt(uint8_t *dst, uint8_t *src, uint8_t key[2])
{
	cipher(BF_ENCRYPT, dst, src, key[0]);
	cipher(BF_ENCRYPT, dst, dst, key[1]);
}


/* Attacks */

const char report[]   = "Key: \"%c%c\"\tOperations: %d\n";
const char notfound[] = "No keys founds in %d operations\n";

typedef struct {
	uint8_t	enc[8]; /* encrypted text, hash table key */
	uint32_t	key;    /* 8-bit encryption key */
} EK;

/* 
 * Meet-in-the-middle attack 
 */

void
encrypt_all_keys(EK items[], uint8_t *plaintext, int size)
{
	int i;

	#pragma omp parallel for private(i)
	for (i = 0; i < size; i++) {
		items[i].key = i;
		cipher(BF_ENCRYPT, items[i].enc, plaintext, i);
	}
}

int compare(const void *i, const void *j)
{
	return memcmp(((EK *)i)->enc,((EK *)j)->enc,8);
}

void
meetinthemiddle(uint8_t *plaintext, uint8_t *ciphertext, int size, int nthreads)
{
	EK *items = (EK *)calloc(size, sizeof(EK));
	int i;

	//printf("encrypt_all_keys\n");
	encrypt_all_keys(items, plaintext, size);

	//printf("qsort\n");
	qsort(items, size, sizeof(int), compare);

	//printf("Parallel for (look for the key)\n");
	#pragma omp parallel for private(i) num_threads(nthreads)
	for (i = 0; i < size; i++) {
		int *found;
		EK lookup;

		//printf("Нить № %d",omp_get_thread_num());
		cipher(BF_DECRYPT, lookup.enc, ciphertext, i);
		if ((found = bsearch(&lookup, items, size, sizeof(int), compare)) != NULL) {
			//printf(report, items[*found].key, i, i+size);
		}
	}
	free(items);

}

#define	MEASURE(x) do {                                              \
						                                       \
	clock_t t = clock(); (x);                                    \
	printf("%.3f\n", (clock()-t)/(double)CLOCKS_PER_SEC);  \
} while (0)

int
main()
{
	uint8_t plaintext[8] = "Dear Bob";
	uint8_t key[2] = "Go";
	uint8_t ciphertext[8];
	uint8_t nthreads, i;
	uint32_t size;
	
	double_encrypt(ciphertext, plaintext, key);

	printf("log size nthreads sec\n");
	for (i=8, size=256; i <= 31; i++, size*=2)
		for (nthreads=1; nthreads<=8; nthreads*=2){
			printf("%d %d %d ",i,size,nthreads);
			MEASURE(meetinthemiddle(plaintext, ciphertext,size,nthreads));
		}
	

	return 0;
}
