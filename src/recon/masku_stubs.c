/*
 * masku_stubs.c -- RECONSTRUCTED (reference build only).
 *
 * core2_decon.c can load a "mask updater" plugin (-umask=...) through
 * Priism's libMaskuPlugin, which drags in the Motif window manager library.
 * The reference build does not support plugins; these stubs make the plugin
 * list empty and every load fail, so -umask reports "plugin not found" and
 * everything else works unchanged.
 *
 * ERDecon reference build, 2026-09-21.
 */
#include "masku_plugin.h"
#include <stddef.h>

void
get_masku_plugin_list(
    const char* config_file, const char* search_path,
    AvailableMaskuPlugin** p_list, int* p_count
)
{
    (void) config_file; (void) search_path;
    *p_list = NULL;
    *p_count = 0;
}

void destroy_masku_plugin_list(AvailableMaskuPlugin* list) { (void) list; }

AvailableMaskuPlugin*
search_masku_plugin_list_by_short_name(
    AvailableMaskuPlugin* list, int count, const char* name
)
{
    (void) list; (void) count; (void) name;
    return NULL;
}

MaskuPlugin load_masku_plugin_from_list(AvailableMaskuPlugin* p) { (void) p; return NULL; }

MaskuPlugin
load_masku_plugin_directly(const char* path, const char* intfc_name)
{
    (void) path; (void) intfc_name;
    return NULL;
}

void unload_masku_plugin(MaskuPlugin plugin) { (void) plugin; }

MaskuPluginAttr* get_masku_plugin_attr(MaskuPlugin plugin) { (void) plugin; return NULL; }
