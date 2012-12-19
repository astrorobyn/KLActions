#ifndef __IO_H__
#define __IO_H__

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>

#include "types.h"
#include "vector.h"
#include "isochroneMap.h"

#define KM_S_TO_KPC_MYR 0.001022729843472133

int read_system_from_text_files(const char *posfile, const char *velfile, body ***bs, size_t *n_read);

#endif // __IO_H__
