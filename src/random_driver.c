#include <random_driver.h>
#include <feature.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

struct _random_plugin_list_t {
  struct _random_plugin_list_t *next;
  random_plugin_t plugin;
};

typedef struct _random_plugin_list_t random_plugin_list_t;

static random_plugin_list_t *random_plugin_list = NULL;
static random_plugin_t *random_plugin = NULL;

/* 2 * uint32_t -> uint64_t */
static uint64_t random_generate64_sim(void) {
  uint64_t upper, lower;

  upper = random_plugin->generate32();
  lower = random_plugin->generate32();
  return (upper << 32) | lower;
}

/* random_driver public variable & API */
double random_gauss_cut = 1.0e35;

double random_get_gcut(void) {
  return random_gauss_cut;
}

double random_set_gcut(double gcut) {
  random_gauss_cut = (gcut < RANDOM_GCUT_MIN) ? RANDOM_GCUT_MIN : gcut;

  return random_gauss_cut;
}

int random_register(const random_plugin_t *plugin) {
  random_plugin_list_t *new, *node, *priv;

  if(plugin == NULL) return -1;

  /* Check random_plugin_t */
  if((RANDOM_ABI_MAJOR(plugin->abi) != RANDOM_ABI_MAJOR(RANDOM_ABI_VERSION)) ||
     (RANDOM_ABI_MINOR(plugin->abi)  < RANDOM_ABI_MINOR(RANDOM_ABI_VERSION))) {
    fprintf(stderr, "random_register(): "
	    "plugin ABI version is mismatched!"
	    "(Required %lu.%lu  Provided %lu.%lu)\n",
	    RANDOM_ABI_MAJOR(RANDOM_ABI_VERSION),
	    RANDOM_ABI_MINOR(RANDOM_ABI_VERSION),
	    RANDOM_ABI_MAJOR(plugin->abi),
	    RANDOM_ABI_MINOR(plugin->abi));
    return -1;
  }

  if(plugin->name == NULL) {
    fprintf(stderr, "random_register(): plugin name is NULL!\n");
    return -1;
  }

  if(plugin->dump == NULL || plugin->restore == NULL) {
    fprintf(stderr, "random_register(): "
	    "plugin[%s] does not support dump/restore function!\n",
	    plugin->name);
    return -1;
  }

  if(plugin->is_supported == NULL
     || plugin->generate == NULL || plugin->generate32 == NULL) {
    fprintf(stderr, "random_register(): "
	    "plugin[%s] does not support generate function family!\n",
	    plugin->name);
    return -1;
  }

  if(plugin->is_supported(RANDOM_UNIFORM0) != 0) {
    fprintf(stderr, "random_register(): "
	    "plugin[%s] does not support uniform distribution [0,1)!\n",
	    plugin->name);
    return -1;
  }

  if(plugin->is_supported(RANDOM_GAUSS_NOCUT) != 0) {
    fprintf(stderr, "random_register(): "
	    "plugin[%s] does not support Gaussian distribution without GCUT!\n",
	    plugin->name);
    return -1;
  }

  if(plugin->is_supported(RANDOM_GAUSS) != 0) {
    fprintf(stderr, "random_register(): "
	    "plugin[%s] does not support Gaussian distribution with GCUT!\n",
	    plugin->name);
    return -1;
  }

  new = malloc(sizeof(random_plugin_list_t));
  if(new == NULL)
    return -1;

  /* Copy plugin description */
  new->plugin.name		= strdup(plugin->name);
  if(new->plugin.name == NULL) {
    free(new);
    return -1;
  }
  if(plugin->sequence != NULL) {
    new->plugin.sequence	= strdup(plugin->sequence);
    if(new->plugin.sequence == NULL) {
      free(new->plugin.name);
      free(new);
      return -1;
    }
  } else {
    new->plugin.sequence	= NULL;
  }

  new->plugin.abi		= plugin->abi;
  new->plugin.version		= plugin->version;
  new->plugin.state_size	= plugin->state_size;
  new->plugin.dump		= plugin->dump;
  new->plugin.restore		= plugin->restore;
  new->plugin.is_supported	= plugin->is_supported;
  new->plugin.generate		= plugin->generate;
  new->plugin.generate32	= plugin->generate32;
  new->plugin.generate64	= plugin->generate64;
  if(new->plugin.generate64 == NULL)
    new->plugin.generate64	= random_generate64_sim;

  /* Check already installed */
  priv = NULL;
  for(node = random_plugin_list; node != NULL; priv = node, node = node->next)
    if(strcmp(node->plugin.name, new->plugin.name) == 0
       && ((   node->plugin.sequence == NULL && new->plugin.sequence == NULL)
	   || (node->plugin.sequence != NULL && new->plugin.sequence != NULL
	       && strcmp(node->plugin.sequence, new->plugin.sequence) == 0))) {
      if(node->plugin.version >= new->plugin.version) {
	fprintf(stderr, "random_register(): "
		"plugin[%s%s%s] version=%" PRIu32 " is already installed"
		"(version=%" PRIu32 " is required to register)\n",
		node->plugin.name,
		node->plugin.sequence == NULL ? "" : "/",
		node->plugin.sequence == NULL ? "" : node->plugin.sequence,
		node->plugin.version, new->plugin.version);
	free(new->plugin.sequence);
	free(new->plugin.name);
	free(new);
	return -1;
      } else {
	fprintf(stderr, "random_register(): "
		"upgrade plugin[%s%s%s] from version=%" PRIu32
		" to %" PRIu32 "\n",
		node->plugin.name,
		node->plugin.sequence == NULL ? "" : "/",
		node->plugin.sequence == NULL ? "" : node->plugin.sequence,
		node->plugin.version, new->plugin.version);
	/* Replace node by new entry */
	new->next = node->next;
	if(priv != NULL)
	  priv->next = new;
	else
	  random_plugin_list = new;

	if(random_plugin == &node->plugin)
	  random_plugin = &new->plugin;

	free(node->plugin.sequence);
	free(node->plugin.name);
	free(node);
	return 0;
      }
    }

  /* Insert new entry to plugin list */
  new->next = random_plugin_list;
  random_plugin_list = new;

  if(random_plugin == NULL)
    random_plugin = &new->plugin;

  return 0;
}

int random_select(const char *name, const char *sequence) {
  random_plugin_list_t *node;

  if(name == NULL) return -1;

  if(random_plugin != NULL && strcmp(random_plugin->name, name) == 0)
    return 0;

  for(node = random_plugin_list; node != NULL; node = node->next)
    if(strcmp(node->plugin.name, name) == 0
       && (sequence == NULL
	   || (node->plugin.sequence != NULL
	       && strcmp(node->plugin.sequence, sequence) == 0))) {
      random_plugin = &node->plugin;
      return 0;
    }

  return -1;
}

const char *random_name(void) {
  if(random_plugin == NULL)
    return NULL;

  return random_plugin->name;
}

const char *random_sequence(void) {
  if(random_plugin == NULL)
    return NULL;

  return random_plugin->sequence;
}

void random_scan_plugin_list(void (*func)(random_plugin_t *, void *),
			     void *param) {
  random_plugin_list_t *node;

  for(node = random_plugin_list; node != NULL; node = node->next)
    func(&(node->plugin), param);
}

size_t random_size(void) {
  random_plugin_list_t *node;
  size_t max;

  for(node = random_plugin_list, max = 0; node != NULL; node = node->next)
    if(node->plugin.state_size > max)
      max = node->plugin.state_size;

  return max;
}

size_t random_dump(uint32_t *buffer, size_t length) {
  if(random_plugin == NULL) return -1;

  return random_plugin->dump(buffer, length);
}

int random_restore(const uint32_t *buffer, size_t length) {
  if(random_plugin == NULL) return -1;

  random_plugin->restore(buffer, length);
  return 0;
}

int    random_check_mode(int mode) {
  if(random_plugin == NULL) return -2;

  return random_plugin->is_supported(mode);
}

int random_generate(int mode, size_t length, double *array) {
  if(random_plugin == NULL) return -1;

  if(length > 0) random_plugin->generate(mode, length, array);
  return 0;
}

uint32_t random_generate32(void) {
  if(random_plugin == NULL) return 0;
  return random_plugin->generate32();
}

uint32_t random_generate64(void) {
  if(random_plugin == NULL) return 0;
  return random_plugin->generate64();
}

int init_framework_Random(void) {
  feature_provide("Random/Framework",
		  FEATURE_VERSION(1, 7), FEATURE_VERSION(1, 6));

  return 0;
}

/* End of File */
