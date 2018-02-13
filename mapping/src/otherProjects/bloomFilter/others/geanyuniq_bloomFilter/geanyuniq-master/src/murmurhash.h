// This file is a part of MapSplice3. Please refer to LICENSE.TXT for the LICENSE
/*
 * Minimalistic header for the MurmurHash2 function. This is in the public domain.
*/

#include <glib.h>

extern guint32	MurmurHash2(gconstpointer key, gsize len, guint32 seed);
