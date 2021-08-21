#ifndef IGNORE_FRAME

// =========================================================================
// Rahmenprogramm fuer Hausarbeit Betriebssysteme SS21
// =========================================================================

#include <time.h>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <pthread.h>

// =========================================================================
// Zeitmessung
// =========================================================================

struct timespec getTimeSpecMonotonic()
{
  struct timespec ts;
  assert(!clock_gettime(CLOCK_MONOTONIC, &ts));
  return ts;
}

// https://www.gnu.org/
//   software/libc/manual/html_node/Calculating-Elapsed-Time.html
// adapted, computes x - y for timespec instead of timeval
int
timespec_subtract(struct timespec xx, 
		  struct timespec yy,
		  struct timespec *result)
{
  struct timespec x = xx, y = yy;
  /* Perform the carry for the later subtraction by updating y. */
  if (x.tv_nsec < y.tv_nsec) {
    int nsec = (y.tv_nsec - x.tv_nsec) / 1000000000L + 1;
    y.tv_nsec -= 1000000000L * nsec;
    y.tv_sec += nsec;
  }
  if (x.tv_nsec - y.tv_nsec > 1000000000L) {
    int nsec = (x.tv_nsec - y.tv_nsec) / 1000000000L;
    y.tv_nsec += 1000000000L * nsec;
    y.tv_sec -= nsec;
  }
  /* Compute the time remaining to wait. tv_nsec is certainly positive. */
  result->tv_sec = x.tv_sec - y.tv_sec;
  result->tv_nsec = x.tv_nsec - y.tv_nsec;
  /* Return 1 if result is negative. */
  return x.tv_sec < y.tv_sec;
}

// convert to us (intended for results of timespec_subtract!)
double
timespec_usec(struct timespec x)
{
  return (double) (1E6 * x.tv_sec + x.tv_nsec / 1E3);
}

double
timeSpecDiffUsec(struct timespec x,
		 struct timespec y)
{
  struct timespec diff;
  timespec_subtract(x, y, &diff);
  return timespec_usec(diff);
}

// =========================================================================
// sequentielles MSB-Radix-Sort
// =========================================================================

// Datentyp fuer Indices in Datenfeld
typedef int64_t SortIndex;
// vorzeichenloser Typ fuer 32-Bit-Key
typedef uint32_t UInt;

// -------------------------------------------------------------------------
// Hilfsfunktionen
// -------------------------------------------------------------------------

void
printFloatData(SortIndex num, UInt *d)
{
  for (SortIndex i = 0; i < num; i++) {
    float f;
    memcpy(&f, d + i, 4);
    printf("%ld\t%g\n", i, f);
  }
}


// -------------------------------------------------------------------------
// bitSorter
// -------------------------------------------------------------------------

// up cond==0 isZero
//  0       0      1
//  0       1      0
//  1       0      0
//  1       1      1
inline int
isZero(UInt up, UInt cond)
{
  return up == (cond == 0);
}

inline void
swap(UInt *a, UInt *b)
{
  UInt tmp = *a;
  *a = *b;
  *b = tmp;
}

// sortiert Datenbereich in d von left bis right (inklusive) in
// Richtung up fuer Bitnummer bitNo
// liefert Start-Index des rechten Bereichs
SortIndex
bitSorter(UInt *d, UInt up, int bitNo, SortIndex left, SortIndex right)
{
  SortIndex l = left, r = right;
  UInt bitMask = ((UInt) 1) << bitNo;
  while (1) {
    // linken Index verschieben
    while ((l <= r) && isZero(up, d[l] & bitMask)) l++;
    // rechten Index verschieben
    while ((l <= r) && !isZero(up, d[r] & bitMask)) r--;
    // Indices ueberkreuzen sich -> Ende
    if (l > r) break;
    // Elemente vertauschen
    swap(d + l, d + r);
  }
  // an diesem Punkt is l = r + 1 (ueberkreuzt rl)
  return l;
}

// -------------------------------------------------------------------------
// Rekursion
// -------------------------------------------------------------------------

void
radixRecursion(UInt *d, UInt up, int bitNo, int lowestBitNo,
	       SortIndex left, SortIndex right)
{
  // printf("radixRecursion %p %u %d %d %ld %ld\n",
  //	 d, up, bitNo, lowestBitNo, left, right);
  if (right <= left) return;
  SortIndex split = bitSorter(d, up, bitNo, left, right);
  // printFloatData(right + 1 - left, d + left);
  bitNo--;
  if (bitNo >= lowestBitNo) {
    radixRecursion(d, up, bitNo, lowestBitNo, left, split - 1);
    radixRecursion(d, up, bitNo, lowestBitNo, split, right);
  }
}

void
radixSortSignAbs(UInt *d, UInt up, int highestBitNo, int lowestBitNo,
		 SortIndex left, SortIndex right)
{
  // printf("radixSortSignAbs %p %u %d %d %ld %ld\n",
  //	 d, up, highestBitNo, lowestBitNo, left, right);
  if (right <= left) return;
  int bitNo = highestBitNo;
  SortIndex split = bitSorter(d, 1 - up, bitNo, left, right);
  // printFloatData(right + 1 - left, d + left);
  bitNo--;
  if (bitNo >= lowestBitNo) {
    radixRecursion(d, 0, bitNo, lowestBitNo, left, split - 1);
    radixRecursion(d, 1, bitNo, lowestBitNo, split, right);
  }
}

// =========================================================================
// paralleles MSB-Radix-Sort
// =========================================================================

#endif // IGNORE_FRAME

// *************************************************************************
// vvvvv Bitte veraendern Sie nur den Bereich ab hier!!! vvvvv
// *************************************************************************

// -------------------------------------------------------------------------
// sort function
// -------------------------------------------------------------------------

typedef struct thread_P {
  unsigned nthreads;
  UInt *d;
  UInt up; 
  int bitNo; 
  int lowestBitNo;
	SortIndex left; 
  SortIndex right;
} rs_param;

void print_rs_param(rs_param *arg, char* message) {
  printf("%s %d %p %u %d %d %ld %ld\n",message,
	  arg->nthreads, arg->d, arg->up, arg->bitNo, arg->lowestBitNo, arg->left, arg->right);
}


void concurrent_radix(unsigned nthreads,UInt *d, UInt up, int bitNo, int lowestBitNo,
	       SortIndex left, SortIndex right);

void* worker_radix(void * initialValues){
  rs_param* parameters = (rs_param*)initialValues;
  concurrent_radix(parameters->nthreads, parameters->d, parameters->up, parameters->bitNo,parameters->lowestBitNo,parameters->left,parameters->right);
  return NULL;
}

void concurrent_radix(unsigned nthreads,UInt *d, UInt up, int bitNo, int lowestBitNo,
	       SortIndex left, SortIndex right) {

  if (right <= left) return;
  SortIndex split = bitSorter(d, up, bitNo, left, right);
  bitNo--;

  if (bitNo >= lowestBitNo) {

    if(nthreads > 1) {

      nthreads/= 2;
      pthread_t thread; 
      rs_param left_param = {nthreads,d,up,bitNo,lowestBitNo,left,split-1};

      pthread_create(&thread, NULL, worker_radix, (void*)&left_param);
      concurrent_radix(nthreads,d,up,bitNo,lowestBitNo,split,right);      
      pthread_join(thread, NULL);
    }
    else {
      concurrent_radix(nthreads,d,up,bitNo,lowestBitNo,left,split-1);
      concurrent_radix(nthreads,d,up,bitNo,lowestBitNo,split,right); 
    }
  }
}

void
parallelRaditxSortSignAbs(unsigned nthreads, UInt *d, UInt up,
			 int highestBitNo, int lowestBitNo,
			 SortIndex left, SortIndex right)
{

  assert(nthreads >= 1);

  if (right <= left) return;
  int bitNo = highestBitNo;
  SortIndex split = bitSorter(d, 1 - up, bitNo, left, right);

  bitNo--;
  if (bitNo >= lowestBitNo) {

    if(nthreads > 1) {

      nthreads/= 2;    
      pthread_t thread;         
      rs_param left_param = {nthreads,d,0,bitNo,lowestBitNo,left,split-1};

      pthread_create(&thread, NULL, worker_radix, (void*)&left_param);
      concurrent_radix(nthreads,d,1,bitNo,lowestBitNo,split,right);      
      pthread_join(thread, NULL);
    }
    else {
      concurrent_radix(nthreads,d,0,bitNo,lowestBitNo,left,split-1);
      concurrent_radix(nthreads,d,1,bitNo,lowestBitNo,split,right); 
    }
  }  
}

// *************************************************************************
// ^^^^^ Bitte veraendern Sie nur den Bereich bis hier!!! ^^^^^
// *************************************************************************

#ifndef IGNORE_FRAME

// =========================================================================
// Anwendung
// =========================================================================

// -------------------------------------------------------------------------
// Erzeugung von Zufallszahlen
// -------------------------------------------------------------------------

void randBytes(uint8_t b[4])
{
  for (int i = 0; i < 4; i++)
    b[i] = rand() & 0xff;
}

float randFloat()
{
  float f;
  uint8_t b[4];
  do {
    randBytes(b);
    memcpy(&f, b, 4);
  } while (!isfinite(f));
  return f;
}

UInt*
generateFloatData(SortIndex num)
{
  UInt *d = malloc(num * sizeof(UInt));
  assert(d);
  for (SortIndex i = 0; i < num; i++) {
    float f = randFloat();
    memcpy(d + i, &f, 4);
  }
  return d;
}

// -------------------------------------------------------------------------
// Test, ob Daten sortiert sind
// -------------------------------------------------------------------------

int
floatsAreSorted(SortIndex num, UInt *d, UInt up)
{
  for (SortIndex i = 1; i < num; i++) {
    float l, r;
    memcpy(&l, d + (i - 1), 4);
    memcpy(&r, d + i, 4);
    if (up ? (l > r) : (l < r))
      return 0;
  }
  return 1;    
}

// =========================================================================
// main
// =========================================================================

int
main(int argc, char *argv[])
{
  assert(sizeof(SortIndex) == 8);
  assert(RAND_MAX > 255);
  // genuegend Parameter vorhanden?
  if (argc != 5) {
    fprintf(stderr, "parameter: <num> <up> <seed> <nthreads>\n");
    exit(-1);
  }
  // Parameter einlesen und ausgeben
  SortIndex num = atol(argv[1]);
  UInt up = atoi(argv[2]);
  unsigned seed = atoi(argv[3]);
  unsigned nthreads = atoi(argv[4]);
  printf("PARAMETER: num = %ld, up = %u, seed = %u, nthreads = %u\n",
	 num, up, seed, nthreads);
  // Zufallszahlen-Generator initialisieren
  srand(seed);
  // Datenfeld erzeugen
  UInt *d = generateFloatData(num);
  // Zeitmessung starten
  struct timespec t0SortMonotonic =
    getTimeSpecMonotonic();
  // Sortierung starten
  if (nthreads == 0) {
    // sequential sorter
    radixSortSignAbs(d, up, 31, 0, 0, num - 1);
  } else {
    // parallel sorter
    parallelRadixSortSignAbs(nthreads, d, up, 31, 0, 0, num - 1);
  }
  // Zeitmessung beenden
  double dtSortMonotonic =
    timeSpecDiffUsec(getTimeSpecMonotonic(), t0SortMonotonic);
  printf("ZEIT: %f us, %f us pro element\n",
	 dtSortMonotonic, dtSortMonotonic / num);
  // Test auf Sortierheit
  // printFloatData(num, d);
  int sorted = floatsAreSorted(num, d, up);
  printf("SORTIERT: %d\n", sorted);
  if (!sorted)
    printf("DATEN SIND NICHT SORTIERT!!!\n");
  return (sorted ? 0 : -1);
}

#endif // IGNORE_FRAME