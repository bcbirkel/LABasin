double  Init_Surface(FILE *fp, surface_t *currentsurface);
double  ShiftDatum_Surface( surface_t *currentsurface,double datum);
int Init_Vel_Profiles(FILE *fp,geologicunit_t *current);
int Init_Vel_Profiles_Binary(FILE *fp,geologicunit_t *current);
int Init_Velocity_Distribution(FILE *fpcontrol,geologicunit_t *current,int id);
void Init_Database(database_t *current, char *databasePath,int visualcheckingoffon);
