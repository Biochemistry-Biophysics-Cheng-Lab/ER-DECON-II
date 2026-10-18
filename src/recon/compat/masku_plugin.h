/*
 * masku_plugin.h -- ERDecon compatibility header.
 *
 * The types and calls of Priism's mask-updater plugin interface that
 * core2_decon.c refers to.  Plugins are not supported in the ERDecon build
 * (recon/masku_stubs.c makes every lookup fail), so this only has to be
 * source compatible; Motif's Widget becomes void*.  ERDecon, 2026-09-23.
 */
#ifndef ERDECON_COMPAT_MASKU_PLUGIN_H
#define ERDECON_COMPAT_MASKU_PLUGIN_H

#include <stddef.h>

struct IVEPluginImpl;
typedef struct IVEPluginImpl* MaskuPlugin;

typedef struct AvailableMaskuPlugin {
    char* name;
    char* short_name;
    char* path;
    char* interface_name;
} AvailableMaskuPlugin;

struct MaskuPluginAttr;

typedef struct MaskuUpdaterClosure {
    struct MaskuPluginAttr* attr;
    void* options;
} MaskuUpdaterClosure;

typedef struct MaskuPluginAttrExt {
    struct MaskuPluginAttrExt* next;
    void* data;
    int type_code;
} MaskuPluginAttrExt;

typedef struct MaskuPluginAttr {
    void* (*options_creator)(struct MaskuPluginAttr* attr);
    void (*options_destroyer)(struct MaskuPluginAttr* attr, void* options);
    int (*arg_parser)(struct MaskuPluginAttr* attr, const char* arg, void* options);
    void (*usage_printer)(struct MaskuPluginAttr* attr,
        void (*printer)(void* completion, const char* fmt, ...), void* completion);
    void (*options_summarizer)(struct MaskuPluginAttr* attr, void* options,
        void (*printer)(void* completion, const char* fmt, ...), void* completion);
    int (*options_menu_creator)(const char* title_prefix, void* parent,
        struct MaskuPluginAttr* attr, void* options);
    int (*thread_count_setter)(struct MaskuPluginAttr* attr, void* options, int count);
    int (*pass_initializer)(struct MaskuPluginAttr* attr, void* options);
    void (*pass_finalizer)(struct MaskuPluginAttr* attr, void* options);
    int (*mask_updater_2d)(int iter, int nx, int ny, const double* guess,
        double* mask, double* work, void* closure);
    int (*mask_updater_3d)(int iter, int nx, int ny, int nz, const double* guess,
        double* mask, double* work, void* closure);
    MaskuPluginAttrExt* extensions;
    void (*error_handler)(void* arg, const char* fmt, ...);
    void* error_handler_arg;
    void* p_internal;
} MaskuPluginAttr;

#ifdef __cplusplus
extern "C" {
#endif
void get_masku_plugin_list(const char* config_file, const char* search_path,
                           AvailableMaskuPlugin** p_list, int* p_count);
void destroy_masku_plugin_list(AvailableMaskuPlugin* list);
AvailableMaskuPlugin* search_masku_plugin_list_by_short_name(
    AvailableMaskuPlugin* list, int count, const char* name);
MaskuPlugin load_masku_plugin_from_list(AvailableMaskuPlugin* p_available);
MaskuPlugin load_masku_plugin_directly(const char* path, const char* intfc_name);
void unload_masku_plugin(MaskuPlugin plugin);
MaskuPluginAttr* get_masku_plugin_attr(MaskuPlugin plugin);
#ifdef __cplusplus
}
#endif

#endif
