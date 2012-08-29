#include "io.h"


int read_system_from_text_files(const char *posfile, const char *velfile, body ***bs, size_t *n_read) {
  FILE *fp, *fv;
    
    //open files for reading 
    fp = fopen(posfile, "r");
    fv = fopen(velfile, "r");
    
    if ((fp != NULL)&&(fv!=NULL)) {
        size_t rstat;
        
        rstat = fscanf(fp, "%ld\n", n_read);
        
        if (rstat != 1) {
            int cstat;
            
            cstat = fclose(fp);
	    fclose(fv);
            if (cstat == EOF) {
                fprintf(stderr, "could not read system size OR close file %s", posfile);
				return 1;
            } else {
                fprintf(stderr, "could not read system size from file %s", posfile);
				return 1;
            }
        } else {
	  size_t i,j;
            
            *bs = (body**) malloc(*n_read * sizeof(body *));
            
            for (i = 0; i < *n_read; i++) {
                body *b;
                size_t rstat;
                
                b = (body*) malloc(sizeof(body));
                
		b->id = i;
		b->m = 1.0;

                rstat = fscanf(fp,"%lg %lg %lg",&(b->q[0]),&(b->q[1]),&(b->q[2]));
		rstat += fscanf(fv,"%lg %lg %lg",&(b->p[0]),&(b->p[1]),&(b->p[2]));
		for(j=0;j<3;j++) b->p[j] *= KM_S_TO_KPC_MYR;

                if (rstat != 6) {
                    int cstat;
                    
                    cstat = fclose(fp);
		    cstat |= fclose(fv);
                    
                    if (cstat == EOF) {
                        fprintf(stderr, "could not read body number %ld OR close files\n", i);
						return 1;
                    } else {
                        fprintf(stderr, "could not read body number %ld from files", i);
						return 1;
                    }
                } else {
                    (*bs)[i] = b;
                }
                
             }
        }
    } else {
        fprintf(stderr, "could not open files for reading");
		return 1;
    }
    
    return 0;
}
