/* Fortran query API for buildinfo framework */
#include <buildinfo.h>
#include <sim/sad_f2c.h>

#include <stdbool.h>
#include <stdio.h>
#include <string.h>

integer buildinfo_get_string_(const_character key, character buf,
			      ftnlen key_len, ftnlen buf_len) {
  const buildinfo_t *entry;

  int key_length = len_trim(key, key_len);
  {
    char key_buf[key_length + 1];
    memcpy(key_buf, key, key_length);
    key_buf[key_length] = '\0';

    entry = buildinfo_search(key_buf);
  }

  if(buildinfo_string_entryQ(entry)) {
    const char *info_string = buildinfo_extract_string(entry);
    if(info_string) {
      size_t info_length = strlen(info_string);
      if(buf_len >= info_length) {
	strncpy(buf, info_string, info_length);
	for(int i = info_length; i < buf_len; i++) {
	  buf[i] = ' ';
	}
	return info_length;
      }
    }
  } else if(buildinfo_integer_entryQ(entry)) {
    char info_string[buf_len + 1];
    int info_length = snprintf(info_string, buf_len + 1,
			       "%d", buildinfo_extract_integer(entry));
    if(buf_len >= info_length) {
      strncpy(buf, info_string, info_length);
      for(int i = info_length; i < buf_len; i++) {
	buf[i] = ' ';
      }
      return info_length;
    }
  }

  return 0;
}

logical buildinfo_get_integer_(const_character key, integer *ret,
			       ftnlen key_len) {
  const buildinfo_t *entry;

  int key_length = len_trim(key, key_len);
  {
    char key_buf[key_length + 1];
    memcpy(key_buf, key, key_length);
    key_buf[key_length] = '\0';

    entry = buildinfo_search(key_buf);
  }

  if(buildinfo_integer_entryQ(entry)) {
    *ret = buildinfo_extract_integer(entry);
    return true;
  }

  return false;
}

/* End of File */
