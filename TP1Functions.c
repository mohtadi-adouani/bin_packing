#include "TP1Functions.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <strings.h>
#include <sys/time.h>
#include <stdio.h>
#include <ilcplex/cplex.h>

int read_TP1_instance(FILE *fin, dataSet *dsptr)
{
	int rval = 0;

	//Taille des boites
	int V;
	//Nombre d'objets
	int n;
	fscanf(fin, "%d,%d\n", &n, &V);
	dsptr->V = V;
	dsptr->n = n;
	dsptr->size = (int *)malloc(sizeof(int) * n);

	int i;
	for (i = 0; i < n; i++)
		fscanf(fin, "%d\n", &(dsptr->size[i]));

	fprintf(stderr, "Instance file read, each bin is %d long and there is %d items of lengths:\n",
			V, n);
	for (i = 0; i < n; i++)
		fprintf(stderr, "%d\t", dsptr->size[i]);
	fprintf(stderr, "\n");

	return rval;
}

int TP1_solve_exact(dataSet *dsptr)
{
	int rval = 0;

	IP_problem *ip_prob_ptr = &(dsptr->master);
	ip_prob_ptr->env = NULL;
	ip_prob_ptr->lp = NULL;
	ip_prob_ptr->env = CPXopenCPLEX(&rval);
	if (rval)
		fprintf(stderr, "ERROR WHILE CALLING CPXopenCPLEX\n");
	if (ip_prob_ptr->env == NULL)
	{
		char errmsg[1024];
		fprintf(stderr, "Could not open CPLEX environment.\n");

		CPXgeterrorstring(ip_prob_ptr->env, rval, errmsg);
		fprintf(stderr, "%s", errmsg);
		exit(0);
	}

	//We create the MIP problem
	ip_prob_ptr->lp = CPXcreateprob(ip_prob_ptr->env, &rval, "TP1");
	if (rval)
		fprintf(stderr, "ERROR WHILE CALLING CPXcreateprob\n");

	rval = CPXsetintparam(ip_prob_ptr->env, CPX_PARAM_DATACHECK, CPX_ON);
	rval = CPXsetintparam(ip_prob_ptr->env, CPX_PARAM_SCRIND, CPX_ON);

	int n = dsptr->n;
	int *size = dsptr->size;
	int V = dsptr->V;
	int nv = n + n * n;

	//We fill our arrays
	//Memory
	ip_prob_ptr->nv = nv;
	ip_prob_ptr->x = (double *)malloc(sizeof(double) * nv);
	ip_prob_ptr->cost = (double *)malloc(sizeof(double) * nv);
	ip_prob_ptr->c_type = (char *)malloc(sizeof(char) * nv);
	ip_prob_ptr->up_bound = (double *)malloc(sizeof(double) * nv);
	ip_prob_ptr->low_bound = (double *)malloc(sizeof(double) * nv);
	ip_prob_ptr->var_name = (char **)malloc(sizeof(char *) * nv);

	int i, j, id = 0;
	//Structures keeping the index of each variable
	int *id_y_i = (int *)malloc(sizeof(int) * n);
	int **id_x_ij = (int **)malloc(sizeof(int *) * n);
	for (i = 0; i < n; i++)
		id_x_ij[i] = (int *)malloc(sizeof(int) * n);

	//First the variables yi (bin #i used or not)
	for (i = 0; i < n; i++)
	{
		//We keep the id
		id_y_i[i] = id;

		//We generate the variable attributes
		ip_prob_ptr->x[id] = 0;
		ip_prob_ptr->cost[id] = 1;
		ip_prob_ptr->c_type[id] = 'B';
		ip_prob_ptr->up_bound[id] = 1;
		ip_prob_ptr->low_bound[id] = 0;
		ip_prob_ptr->var_name[id] = (char *)malloc(sizeof(char) * 1024);
		snprintf(ip_prob_ptr->var_name[id],
				 1024,
				 "y_i%d",
				 i);
		id++;
	}

	// Initialisation pour chaque Xij
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			//We keep the id
			id_x_ij[i][j] = id;

			//We generate the variable attributes
			ip_prob_ptr->x[id] = 0;		   // On lui done la valeur 0 (valeur arbitraire)
			ip_prob_ptr->cost[id] = 0;	   // coeff dans la fonction objective
			ip_prob_ptr->c_type[id] = 'B'; //Type boolean
			//--------Domaine de definition-----
			ip_prob_ptr->up_bound[id] = 1;
			ip_prob_ptr->low_bound[id] = 0;
			//----------------------------------
			//Nom de chaque variable
			ip_prob_ptr->var_name[id] = (char *)malloc(sizeof(char) * 1024);
			snprintf(ip_prob_ptr->var_name[id],
					 1024,
					 "x_i%d_j%d",
					 i, j);

			id++;
		}
	}

	/********************************************/
	/******* 	     VARIABLES x_ij 	  *******/
	/********************************************/

	//Ajout d'une variable initialisee a Cplex
	rval = CPXnewcols(ip_prob_ptr->env, ip_prob_ptr->lp,
					  nv,
					  ip_prob_ptr->cost,
					  ip_prob_ptr->low_bound,
					  ip_prob_ptr->up_bound,
					  ip_prob_ptr->c_type,
					  ip_prob_ptr->var_name);
	if (rval)
		fprintf(stderr, "CPXnewcols returned errcode %d\n", rval);

	//Allocation de memoire pour les outils
	//Constraints part
	ip_prob_ptr->rhs = (double *)malloc(sizeof(double));
	ip_prob_ptr->sense = (char *)malloc(sizeof(char));
	ip_prob_ptr->rmatbeg = (int *)malloc(sizeof(int));
	ip_prob_ptr->nz = n + 1;
	ip_prob_ptr->rmatind = (int *)malloc(sizeof(int) * nv);
	ip_prob_ptr->rmatval = (double *)malloc(sizeof(double) * nv);
	ip_prob_ptr->const_name = (char **)malloc(sizeof(char *));
	ip_prob_ptr->const_name[0] = (char *)malloc(sizeof(char) * 1024);

	ip_prob_ptr->rmatbeg[0] = 0;
	ip_prob_ptr->rhs[0] = 0;
	ip_prob_ptr->sense[0] = 'L'; //Inferieur ou egale

	//Partie gauche des contraintes
	for (i = 0; i < n; i++)
	{

		snprintf(ip_prob_ptr->const_name[0], 1024, "bin_sac_i%d", i);
		id = 0;
		ip_prob_ptr->rmatind[id] = id_y_i[i];
		ip_prob_ptr->rmatval[id] = -V;
		id++;

		for (j = 0; j < n; j++)
		{
			ip_prob_ptr->rmatind[id] = id_x_ij[i][j];
			ip_prob_ptr->rmatval[id] = size[j];
			id++;
		}

		//Ajout de la contraintes a Cplex
		rval = CPXaddrows(ip_prob_ptr->env, ip_prob_ptr->lp,
						  0,
						  1,
						  n + 1,
						  ip_prob_ptr->rhs,
						  ip_prob_ptr->sense,
						  ip_prob_ptr->rmatbeg,
						  ip_prob_ptr->rmatind,
						  ip_prob_ptr->rmatval,
						  NULL,
						  ip_prob_ptr->const_name);
		if (rval)
			fprintf(stderr, "CPXaddrows 1 returned errcode %d\n", rval);
	}

	ip_prob_ptr->rhs[0] = 1;

	ip_prob_ptr->sense[0] = 'E'; //egale
	ip_prob_ptr->rmatbeg[0] = 0;

	ip_prob_ptr->rmatval[n] = 0;
	ip_prob_ptr->rmatind[n] = 0;

	for (i = 0; i < n; i++)
	{

		snprintf(ip_prob_ptr->const_name[0], 1024, "objectUsed_i%d", i);
		for (j = 0; j < n; j++)
		{
			ip_prob_ptr->rmatval[j] = 1;
			ip_prob_ptr->rmatind[j] = id_x_ij[j][i];
		}
		//Ajout de la contraintes a Cplex
		rval = CPXaddrows(ip_prob_ptr->env, ip_prob_ptr->lp,
						  0,
						  1,
						  n + 1,
						  ip_prob_ptr->rhs,
						  ip_prob_ptr->sense,
						  ip_prob_ptr->rmatbeg,
						  ip_prob_ptr->rmatind,
						  ip_prob_ptr->rmatval,
						  NULL,
						  ip_prob_ptr->const_name);
		if (rval)
		{
			fprintf(stderr, "CPXaddrows 2 returned errcode %d\n", rval);
		}
	}

	//We write the problem for debugging purposes, can be commented afterwards
	rval = CPXwriteprob(ip_prob_ptr->env, ip_prob_ptr->lp, "bin_packing.lp", NULL);
	if (rval)
		fprintf(stderr, "CPXwriteprob returned errcode %d\n", rval);

	//We solve the model
	rval = CPXmipopt(ip_prob_ptr->env, ip_prob_ptr->lp);
	if (rval)
		fprintf(stderr, "CPXmipopt returned errcode %d\n", rval);

	rval = CPXsolwrite(ip_prob_ptr->env, ip_prob_ptr->lp, "bin_packing.sol");
	if (rval)
		fprintf(stderr, "CPXsolwrite returned errcode %d\n", rval);

	//We get the objective value
	rval = CPXgetobjval(ip_prob_ptr->env, ip_prob_ptr->lp, &(ip_prob_ptr->objval));
	if (rval)
		fprintf(stderr, "CPXgetobjval returned errcode %d\n", rval);

	//We get the best solution found
	rval = CPXgetobjval(ip_prob_ptr->env, ip_prob_ptr->lp, &(ip_prob_ptr->objval));
	rval = CPXgetx(ip_prob_ptr->env, ip_prob_ptr->lp, ip_prob_ptr->x, 0, nv - 1);
	if (rval)
		fprintf(stderr, "CPXgetx returned errcode %d\n", rval);

	//We display the solution
	double tolerance = 0.0001;
	int remaining;
	for (i = 0; i < n; i++)
	{
		id = id_y_i[i];
		if (ip_prob_ptr->x[id] <= 1 - tolerance)
			continue;
		remaining = V;
		fprintf(stderr, "Bin #%d is used with volume %d and contains:\n", i, V);
		for (j = 0; j < n; j++)
		{
			id = id_x_ij[i][j];
			if (ip_prob_ptr->x[id] <= 1 - tolerance)
				continue;
			remaining -= size[j];
			fprintf(stderr, "\tItem #%d of volume %d (remaining: %d)\n", j, size[j], remaining);
		}
	}

	return rval;
}

//Afficher integer tab
void show_int_tab(int *tab, int n)
{
	for (int i = 0; i < n; i++)
	{
		printf("||Boite %d => place restante %d.\n",i, tab[i]);
	}
	printf("-----------------\n");
}

//Trie un tab entier decroissant
int * sorte_tab_dec(int tab[], int n){
    int temp;
    int size = n * sizeof(int);
    int * temp_tab = (int*)malloc(size);
    memcpy(temp_tab, tab, size);
    for(int i = 0; i < n-1 ; i++ ){
        for(int j = i+1 ; j < n ;j++){
            if (temp_tab[i] < temp_tab[j]){
                temp = temp_tab[i];
                temp_tab[i] = temp_tab[j];
                temp_tab[j] = temp;
            }
        }
    }
    return temp_tab;
}


//FIRST FIT ALGO
int heuriFirstFit(int tab[], int box[], int n)
{
	int rval = 0; //nb de box used 
	int size = n * sizeof(int); //memoire size
	int * temp_box = (int *)malloc(size);
	memcpy(temp_box, box, size); // copie le tableau
	int init = box[0];
	
	// objects iteration
	for (int i = 0; i < n; i++)
	{
		// box iteration
		for (int j = 0; j < n; j++)
		{
			if (tab[i] <= temp_box[j]) // si assez de place
			{	
				printf("Objet %d mis dans la boite %d de taille %d.\n ",i,j,tab[i]);
				temp_box[j] = temp_box[j]  - tab[i]; //On met dans la boite
				break;
			}
		}
	}
	while (temp_box[rval] != init && rval < n)//tant que la derniere boite n'est pas pleine et que nous n'avons pas uttilisee tt les boites.
	{
		rval++;
	}

    fprintf(stderr, "[Heurestiques FIRST FIT] ,Il faut %d boites.\n", rval);
	show_int_tab(temp_box, n);
	
	return rval;
}



//NEXT FIT ALOG
int heuriNextFit(int tab[], int box[], int n)
{
	int rval = 1;
	int temp = 0;
	int size = n * sizeof(int);
	int * temp_box = (int *)malloc(size);
	memcpy(temp_box, box, size);
	// box iteration
	for (int i = 0; i < n; i++)
	{
		// object iteration
		for (int j = temp; j < n; j++)
		{
			if (temp_box[i] >= tab[j]) // si il reste de la place
			{
				temp_box[i] = temp_box[i] - tab[j]; //on met dand la boite (on considere) 
				printf("Objet %d mis dans la boite %d de taille %d.\n ",j,i,tab[i]);
				temp++;
			}
			else // sinon nouvelle boite
			{
				rval++;
				break;
			}
		}
	}
	fprintf(stderr, "[HEURESTIQUE NEXT FIT],Il faut %d boites.\n", rval);
	show_int_tab(temp_box, n);
	
	return rval;
}

int heuriFirstFitDecreasing(int tab[], int box[], int n){
	int rval = 0;

	int size = n * sizeof(int);
	int * temp_box = (int *)malloc(size);
	memcpy(temp_box, box, size);


    // Pour chaque objet
    for (int i = 0; i < n; i++) {
        //Chercher la boite 
		int j;
        for (j = 0; j < rval; j++) {
            if (temp_box[j] >= tab[i]) {//si la boite peut le contenir
                temp_box[j] = temp_box[j] - tab[i]; //le mettre dans la boite
               	printf("Objet %d mis dans la boite %d de taille %d.\n ",i,j,tab[i]);
                break;
            }
        }
 
        //Nouvelle Boite
        if (j == rval) {
            temp_box[rval] = temp_box[i] - tab[i];
            rval++;
        }
       
    }

	printf("[HEURESTIQUE FIRST FIT DECREASING],Il faut %d boites.\n", rval);

	show_int_tab(temp_box,n);
	return rval;
}


int TP1_solve_heuristic(dataSet *dsptr)
{
    int n = dsptr->n;
    
	int * box_v = (int *)malloc(n * sizeof(int));
	for (int i = 0 ; i < n ; i++)
	{
		box_v[i] = dsptr->V;
	}

	//Pour recuperer la valeur dans une variaable il suffit de mettre le resuatat de des fonction dans une variable
	
	printf("################# First Fit #########################\n");
	heuriFirstFit(dsptr->size, box_v, dsptr->n);


	printf("################# Next Fit #########################\n");
	heuriNextFit(dsptr->size, box_v, dsptr->n);


	printf("################# First Fit Decreasing #########################\n");
	int * tabSorted = sorte_tab_dec(dsptr->size, dsptr->n);
	heuriFirstFitDecreasing(tabSorted, box_v, dsptr->n);

	return 0;
}
