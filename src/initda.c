#include <sim/sad_f2c.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

static void wait_prompt(char*);
static int read_line(int, char*);

static int pfp1[2],pfp2[2];
static pid_t pid;

void initda_(integer4 *nord, integer4 *lvec) {
   char str[256];
   char *s1;
   int iprmpt;

   pipe(pfp1);
   pipe(pfp2);
   switch(pid = fork()) {
   case -1:
     printf("fork err\n");
     break;

   case 0:
     printf("forked process\n");
     dup2(pfp1[0], STDIN_FILENO);
     dup2(pfp2[1], STDOUT_FILENO);
     execlp("/proj/oldsad/sad+", NULL);
     printf("chiled: exec error.\n");
     break;

   default:
     wait_prompt(str);
     sprintf(str, "DAinit[%d]\n", *nord);
     write(pfp1[1], str, strlen(str));
     iprmpt = 0;
     do {
       iprmpt = read_line(pfp2[0], str);
       if(!iprmpt) {
	 printf("%s", str);
	 if((s1 = strstr(str, "Length of DA vector =")) != NULL) 
	   *lvec = atoi(s1 + 21);
       }
       else
	 printf("***/DA:***>");
     } while(strstr(str, "%>") == NULL);
     break;
   }
}

void execda_(integer4 *nord) {
  FILE *fpin;
  char str[256];

  fpin = popen("/proj/oldsad/bad", "r");
  while(fgets(str, 256, fpin) != NULL) {
    printf("%s", str);
  }
  pclose(fpin);
}

void DAexit(void) {
  if(pid == 0) return;
  write(pfp1[1], "Exit[]\n", 7);
}

static void wait_prompt(char *str) {
  int iprmpt = 0;
  do {
    iprmpt = read_line(pfp2[0], str);
    if(!iprmpt)
      printf("%s", str);
    else
      printf("***/DA:***>");
  } while(strstr(str, "%>") == NULL);
}

static int read_line(int fp, char *str) {
   int i = 0, iprmpt = 0, ipsearch = 0;
   char c;

   do {
     read(fp, &c, 1);
     str[i] = c; i++;
     if(c == '%') ipsearch = 1;
     if(c == '>' && ipsearch) {
       c = '\n';
       str[i] = c; i++;
       iprmpt=1;
     }
   } while(c != '\n');
   str[i] = '\0';
   return iprmpt;
}

/* End of File */
