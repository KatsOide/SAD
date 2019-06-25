#ifndef CALDATE_H
#define CALDATE_H

struct caldate {
  long year;
  int month;
  int day;
} ;

extern unsigned int caldate_fmt(char *, struct caldate *);
extern unsigned int caldate_scan(char *, struct caldate *);

extern void caldate_frommjd(struct caldate *, long, int *, int *);
extern long caldate_mjd(struct caldate *);
extern void caldate_normalize(struct caldate *);

extern void caldate_easter(struct caldate *);

#endif
