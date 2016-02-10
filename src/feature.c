#include <feature.h>

#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

struct _feature_entry;

struct _feature_entry {
  struct _feature_entry *next;
  uint32_t version;
  uint32_t version_compatible;
  const char *name;
};

typedef struct _feature_entry feature_entry_t;

static feature_entry_t *feature_db = NULL;

bool feature_provide(const char *name,
		     uint32_t version, uint32_t version_compatible) {
  feature_entry_t *entry, *priv, *new;

  if(version_compatible > version) {
    fprintf(stderr, "feature_provide() error: version=%" PRIu32
	    " < version_compatible=%" PRIu32 "\n",
	    version, version_compatible);
    return false;
  }

  priv = NULL;
  for(entry = feature_db; entry != NULL; priv = entry, entry = entry->next)
    if(strcmp(name, entry->name) == 0) {
      if(entry->version > version
	 || (entry->version == version
	     && version_compatible >= entry->version_compatible)) {
	fprintf(stderr, "feature[%s] is already provided as version=%" PRIu32
		"(feature_provide is called with version=%" PRIu32 ")\n",
		name, entry->version, version);
	return false;
      }
      if(version > entry->version) {
	fprintf(stderr, "feature[%s] version is redefined from %" PRIu32
		" to %" PRIu32 ".\n",
		name, entry->version, version);
      }
      if(version_compatible > entry->version_compatible) {
	fprintf(stderr, "feature[%s] version compatibility is limited"
		" from %" PRIu32
		" to %" PRIu32 ".\n",
		name, entry->version_compatible, version_compatible);
      }
      /* Override current entry */
      entry->version = version;
      entry->version_compatible = version_compatible;
      return true;
    }

  new = malloc(sizeof(feature_entry_t));
  if(new == NULL) {
    fprintf(stderr, "Can't allocate memory for new feature_db entry!\n");
    return false;
  }

  new->next = NULL;
  new->version = version;
  new->version_compatible = version_compatible;
  new->name = strdup(name);
  if(new->name == NULL) {
    fprintf(stderr, "Can't allocate memory for new feature_db entry!\n");
    free(new);
    return false;
  }

  if(priv != NULL)
    priv->next = new;
  else
    feature_db = new;

  return true;
}

bool feature_require(const char *name, uint32_t version) {
  feature_entry_t *entry;

  for(entry = feature_db; entry != NULL; entry = entry->next)
    if(strcmp(name, entry->name) == 0)
      return entry->version >= version && version >= entry->version_compatible;

  return false;
}

uint32_t feature_version(const char *name) {
  feature_entry_t *entry;

  for(entry = feature_db; entry != NULL; entry = entry->next)
    if(strcmp(name, entry->name) == 0)
      return entry->version;

  return 0UL;
}

int feature_scan(void (*func)(const char*)) {
  feature_entry_t *entry;
  int count = 0;

  for(entry = feature_db; entry != NULL; entry = entry->next, count += 1)
    func(entry->name);

  return count;
}

int init_framework_Feature(void) {
  feature_provide("Feature/Framework",
		  FEATURE_VERSION(1, 0), FEATURE_VERSION(1, 0));

  return 0;
}

/* End of File */
