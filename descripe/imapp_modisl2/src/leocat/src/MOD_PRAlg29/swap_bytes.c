# include <stdlib.h>

int swapbyte4(void *, int);
void swbyt4_(void *, int *);
void fbyte4_(void *, int *);

/* ===== swapbyte4 ================================================== */
/*
*| Name:
*|    swapbyte4 - Flips 4 byte words into network byte order.
*|
*| Interface:
*|    int
*|    swapbyte4 (void *buf, int n_words)
*|
*| Input:
*|    buf      - Buffer to be flipped.
*|    n_words  - Number of words to flip.
*|
*| Input and Output:
*|    none
*|
*| Output:
*|    none
*|
*| Return values:
*|    0        - Success.
*|   -1        - buf pointer is NULL.
*|
*| Remarks:
*|    Network byte order is big endian.
*|
*| Categories:
*|    converter
*|    utility
*/

int
  swapbyte4 (void *buf, int n_words)
{
  int          *b_ptr = (int *) buf;
  int           nw    = (int) n_words;

  if (b_ptr == (int *) NULL)
  {
    return (-1);
  }

  if (n_words > 0)
  {
    swbyt4_ (b_ptr, &nw);
  }
  return (0);
}

/*
*$ Name:
*$      swbyt4   - Re-orders bytes if internal representation is not
*$                 big-endian.
*$
*$ Interface:
*$      subroutine
*$      swbyt4(integer buf(*), integer n)
*$
*$ Input:
*$      n        - Number of 4 byte swaps to be made if order is not
*$                 big-endian.
*$
*$ Input and Output:
*$      buf      - Array containing bytes to be swapped.
*$
*$ Output:
*$      none
*$
*$ Return values:
*$      none
*$
*$ Remarks:
*$      none
*$
*$ Categories:
*$      utility
*/

void
swbyt4_(void *buf, int *n)
{
        /* determine byte order from first principles */
        union
        {
                char            bytes[sizeof(int)];
                int             word;
        }               q;

        q.word = 1;
        if (q.bytes[3] == 1)
                return;

        /* swap bytes */
        fbyte4_(buf, n);
}

   /* buffer - INTEGER*4 array to switch bytes  */
   /* n      - INTEGER number of 4 byte switches */

void
fbyte4_(void *buffer, int *num)
{
  char *cbuf = (char *)buffer;
  int i, n;

  n = *num;
  for (i = 0; i < n; i++)
  {
    char b;

    b = cbuf[0];
    cbuf[0] = cbuf[3];
    cbuf[3] = b;
    b = cbuf[1];
    cbuf[1] = cbuf[2];
    cbuf[2] = b;
    cbuf += 4;
  }
}

