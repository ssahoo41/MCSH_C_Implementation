// void MaxwellCartesianSphericalHarmonics(const double *x, const double *y, const double *z, const int l, const char *n, const double rCutoff, double *result, const int size)
// {
// 	double *r = malloc( size * sizeof(double));
// 	double *x_hat = malloc( size * sizeof(double));
// 	double *y_hat = malloc( size * sizeof(double));
// 	double *z_hat = malloc( size * sizeof(double));


// 	getRArray(x, y, z, r, size);
// 	divideVector(x, r, x_hat, size);
// 	divideVector(y, r, y_hat, size);
// 	divideVector(z, r, z_hat, size);

// 	double *uncutResult = malloc( size * sizeof(double));
// 	// double *result = malloc( size * sizeof(double));

// 	int i;
// 	// printf("\n============\n");
// 	// for (i = 0; i < size; i++)
// 	// {
// 	// 	printf("x: %10f \t y: %10f \t z: %10f \t x_hat: %10f \t y_hat: %10f \t z_hat: %10f \t r: %10f\n", x[i],y[i],z[i],x_hat[i],y_hat[i],z_hat[i], r[i]);
// 	// }

// 	// for (i = 0; i < size; i++)
// 	// {
// 	// 	double r_calc = sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
// 	// 	printf("x: %10f \t y: %10f \t z: %10f \t r: %10f \t r_calc: %10f \n", x[i],y[i],z[i], r[i],  r_calc);
// 	// }
// 	switch (l) 
// 	{
// 		case 0:
// 			//uncutResult = 1;
			
// 			for ( i = 0; i < size; i++);
// 				uncutResult[i] = 1.0;
// 			break;

// 		case 1:
// 			if (strcmp(n, "100") == 0) 
// 			{
// 				for ( i = 0; i < size; i++)
// 				{
// 					uncutResult[i] = x_hat[i];

// 				}
// 			} 
// 			else if (strcmp(n, "010") == 0)
// 			{
// 				for ( i = 0; i < size; i++)
// 				{
// 					uncutResult[i] = y_hat[i];
// 				}
// 			}
// 			else if (strcmp(n, "001") == 0)
// 			{
// 				for ( i = 0; i < size; i++)
// 				{
// 					uncutResult[i] = z_hat[i];
// 				}
// 			}
// 			else
// 			{
// 				printf("\nWARNING: n is not valid %s \n", n);
// 			}
// 			break;

// 		case 2:
// 			if (strcmp(n, "200") == 0) 
// 			{
// 				// result = 3.0 * x_hat * x_hat - 1.0;
// 				double *temp = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 2, 0, 0, 3.0, temp, size);
// 				addScalarVector(temp, -1.0, uncutResult, size);
// 				free(temp);
// 			} 
// 			else if (strcmp(n, "020") == 0)
// 			{
// 				// result = 3.0 * y_hat * y_hat - 1.0;
// 				double *temp = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 2, 0, 3.0, temp, size);
// 				addScalarVector(temp, -1.0, uncutResult, size);
// 				free(temp);
// 			}
// 			else if (strcmp(n, "002") == 0)
// 			{
// 				// result = 3.0 * z_hat * z_hat - 1.0;
// 				double *temp = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 2, 3.0, temp, size);
// 				addScalarVector(temp, -1.0, uncutResult, size);
// 				free(temp);
// 			}
// 			else if (strcmp(n, "110") == 0)
// 			{
// 				// result = 3.0 * x_hat * y_hat;
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 1, 0, 3.0, uncutResult, size);
// 			}
// 			else if (strcmp(n, "101") == 0)
// 			{
// 				// result = 3.0 * x_hat * z_hat;
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 1, 3.0, uncutResult, size);
// 			}
// 			else if (strcmp(n, "011") == 0)
// 			{
// 				// result = 3.0 * y_hat * z_hat;
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 1, 3.0, uncutResult, size);
// 			}
// 			else
// 			{
// 				printf("\nWARNING: n is not valid %s \n", n);
// 			}
// 			break;

// 		case 3:
// 			if (strcmp(n, "300") == 0) 
// 			{
// 				// result = 15.0 * x_hat * x_hat * x_hat - 9.0 * x_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 3, 0, 0, 15.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 0, 9.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			} 
// 			else if (strcmp(n, "030") == 0)
// 			{
// 				// result = 15.0 * y_hat * y_hat * y_hat - 9.0 * y_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 3, 0, 15.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 0, 9.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "003") == 0)
// 			{
// 				// result = 15.0 * z_hat * z_hat * z_hat - 9.0 * z_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 3, 0, 15.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 0, 9.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "210") == 0)
// 			{
// 				// result = 15.0 * x_hat * x_hat * y_hat - 3.0 * y_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 2, 1, 0, 15.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 0, 3.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "120") == 0)
// 			{
// 				// result = 15.0 * x_hat * y_hat * y_hat - 3.0 * x_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 2, 0, 15.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 0, 3.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "201") == 0)
// 			{
// 				// result = 15.0 * x_hat * x_hat * z_hat - 3.0 * z_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 2, 0, 1, 15.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 1, 3.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "102") == 0)
// 			{
// 				// result = 15.0 * x_hat * z_hat * z_hat - 3.0 * x_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 2, 15.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 0, 3.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "021") == 0)
// 			{
// 				// result = 15.0 * y_hat * y_hat * z_hat - 3.0 * z_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 2, 1, 15.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 1, 3.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "012") == 0)
// 			{
// 				// result = 15.0 * y_hat * z_hat * z_hat - 3.0 * y_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 2, 15.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 0, 3.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "111") == 0)
// 			{
// 				// result = 15.0 * x_hat * y_hat * z_hat;
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 1, 1, 15.0, uncutResult, size);
// 			}
// 			else
// 			{
// 				printf("\nWARNING: n is not valid %s \n", n);
// 			}
// 			break;

// 		case 4:
// 			if (strcmp(n, "400") == 0) 
// 			{
// 				// result = 105.0 * x_hat * x_hat * x_hat * x_hat - 90.0 * x_hat * x_hat + 9.0;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				double *temp3 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 4, 0, 0, 105.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 2, 0, 0, 90.0, temp2, size);
// 				subtractVector(temp1, temp2, temp3, size);
// 				addScalarVector(temp3, 9.0, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 				free(temp3);
// 			} 
// 			else if (strcmp(n, "040") == 0)
// 			{
// 				// result = 105.0 * y_hat * y_hat * y_hat * y_hat - 90.0 * y_hat * y_hat + 9.0;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				double *temp3 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 4, 0, 105.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 2, 0, 90.0, temp2, size);
// 				subtractVector(temp1, temp2, temp3, size);
// 				addScalarVector(temp3, 9.0, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 				free(temp3);
// 			}
// 			else if (strcmp(n, "004") == 0)
// 			{
// 				// result = 105.0 * z_hat * z_hat * z_hat * z_hat - 90.0 * z_hat * z_hat + 9.0;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				double *temp3 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 4, 105.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 2, 90.0, temp2, size);
// 				subtractVector(temp1, temp2, temp3, size);
// 				addScalarVector(temp3, 9.0, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 				free(temp3);
// 			}
// 			else if (strcmp(n, "310") == 0)
// 			{
// 				// result = 105.0 * x_hat * x_hat * x_hat * y_hat - 45.0 * x_hat * y_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 3, 1, 0, 105.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 1, 0, 45.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "130") == 0)
// 			{
// 				// result = 105.0 * x_hat * y_hat * y_hat * y_hat - 45.0 * x_hat * y_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 3, 0, 105.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 1, 0, 45.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "301") == 0)
// 			{
// 				// result = 105.0 * x_hat * x_hat * x_hat * z_hat - 45.0 * x_hat * z_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 3, 0, 1, 105.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 1, 45.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "103") == 0)
// 			{
// 				// result = 105.0 * x_hat * z_hat * z_hat * z_hat - 45.0 * x_hat * z_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 3, 105.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 1, 45.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "031") == 0)
// 			{
// 				// result = 105.0 * y_hat * y_hat * y_hat * z_hat - 45.0 * y_hat * z_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 3, 1, 105.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 1, 45.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "013") == 0)
// 			{
// 				// result = 105.0 * y_hat * z_hat * z_hat * z_hat - 45.0 * y_hat * z_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 3, 105.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 1, 45.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "220") == 0)
// 			{
// 				// result = 105.0 * x_hat * x_hat * y_hat * y_hat - 15.0 * x_hat * x_hat - 15.0 * y_hat * y_hat + 3.0;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				double *temp3 = malloc( size * sizeof(double));
// 				double *temp4 = malloc( size * sizeof(double));
// 				double *temp5 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 2, 2, 0, 105.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 2, 0, 0, 15.0, temp2, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 2, 0, 15.0, temp3, size);
// 				subtractVector(temp1, temp2, temp4, size);
// 				subtractVector(temp4, temp3, temp5, size);
// 				addScalarVector(temp5, 3.0, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 				free(temp3);
// 				free(temp4);
// 				free(temp5);
// 			}
// 			else if (strcmp(n, "202") == 0)
// 			{
// 				// result = 105.0 * x_hat * x_hat * z_hat * z_hat - 15.0 * x_hat * x_hat - 15.0 * z_hat * z_hat + 3.0;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				double *temp3 = malloc( size * sizeof(double));
// 				double *temp4 = malloc( size * sizeof(double));
// 				double *temp5 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 2, 0, 2, 105.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 2, 0, 0, 15.0, temp2, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 2, 15.0, temp3, size);
// 				subtractVector(temp1, temp2, temp4, size);
// 				subtractVector(temp4, temp3, temp5, size);
// 				addScalarVector(temp5, 3.0, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 				free(temp3);
// 				free(temp4);
// 				free(temp5);
// 			}
// 			else if (strcmp(n, "022") == 0)
// 			{
// 				// result = 105.0 * y_hat * y_hat * z_hat * z_hat - 15.0 * y_hat * y_hat - 15.0 * z_hat * z_hat + 3.0;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				double *temp3 = malloc( size * sizeof(double));
// 				double *temp4 = malloc( size * sizeof(double));
// 				double *temp5 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 2, 2, 105.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 2, 0, 15.0, temp2, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 0, 2, 15.0, temp3, size);
// 				subtractVector(temp1, temp2, temp4, size);
// 				subtractVector(temp4, temp3, temp5, size);
// 				addScalarVector(temp5, 3.0, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 				free(temp3);
// 				free(temp4);
// 				free(temp5);
// 			}
// 			else if (strcmp(n, "211") == 0)
// 			{
// 				// result = 105.0 * x_hat * x_hat * y_hat * z_hat - 15.0 * y_hat * z_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 2, 1, 1, 105.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 0, 1, 1, 15.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "121") == 0)
// 			{
// 				// result = 105.0 * x_hat * y_hat * y_hat * z_hat - 15.0 * x_hat * z_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 2, 1, 105.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 0, 1, 15.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else if (strcmp(n, "112") == 0)
// 			{
// 				// result = 105.0 * x_hat * y_hat * z_hat * z_hat - 15.0 * x_hat * y_hat;
// 				double *temp1 = malloc( size * sizeof(double));
// 				double *temp2 = malloc( size * sizeof(double));
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 1, 2, 105.0, temp1, size);
// 				polyXYZArray(x_hat, y_hat, z_hat, 1, 1, 0, 15.0, temp2, size);
// 				subtractVector(temp1, temp2, uncutResult, size);
// 				free(temp1);
// 				free(temp2);
// 			}
// 			else
// 			{
// 				printf("\nWARNING: n is not valid %s \n", n);
// 			}
// 			break;

// 		default:
// 			printf("\nWARNING: l is not valid %d \n", l);
// 	}

// 	//int i;
// 	for (i = 0; i < size; i++)
// 	{
// 		if (r[i] < rCutoff)
// 		{
// 			result[i] = uncutResult[i];
// 		}
// 		else
// 		{
// 			result[i] = 0;
// 		}
// 		//printf("%10f \t %10f \n", uncutResult[i], result[i]);
// 	}


// 	free(uncutResult);
// 	free(r);
// 	free(x_hat);
// 	free(y_hat);
// 	free(z_hat);

// 	// return result;
// }













/*double *refX = malloc( pixelEvalArrSize * sizeof(double));
	double *refY = malloc( pixelEvalArrSize * sizeof(double));
	double *refZ = malloc( pixelEvalArrSize * sizeof(double));

	getCentralCoords(hx, hy, hz, accuracy, refX, refY, refZ);

	int centerX = (stencilDimX - 1)/2;
    int centerY = (stencilDimY - 1)/2;
    int centerZ = (stencilDimZ - 1)/2;

    double *stencil = malloc( stencilDimX * stencilDimY * stencilDimZ * sizeof(double));

	double *tempXArr = malloc( pixelEvalArrSize * sizeof(double));
	double *tempYArr = malloc( pixelEvalArrSize * sizeof(double));
	double *tempZArr = malloc( pixelEvalArrSize * sizeof(double));
	double *tempMCSHResult = malloc( pixelEvalArrSize * sizeof(double));
	double xOffset, yOffset, zOffset;
	int i, j, k, index = 0;
	for (k = 0; k < stencilDimZ; k++){
		for ( j = 0; j < stencilDimY; j++) {
			for ( i = 0; i < stencilDimX; i++) {
				// printf("start %d %d %d\n", i,j,k);
				xOffset = (i-centerX) * hx;
				yOffset = (j-centerY) * hy;
				zOffset = (k-centerZ) * hz;
				//index = k * stencilDimX * stencilDimY + j * stencilDimX + i;

				addScalarVector(refX, xOffset, tempXArr, pixelEvalArrSize);
				addScalarVector(refY, yOffset, tempYArr, pixelEvalArrSize);
				addScalarVector(refZ, zOffset, tempZArr, pixelEvalArrSize);

				applyU2(tempXArr, tempYArr, tempZArr, U, pixelEvalArrSize);

				MaxwellCartesianSphericalHarmonics(tempXArr, tempYArr, tempZArr, l, n, rCutoff, tempMCSHResult, pixelEvalArrSize);

				//stencil[index] = cblas_dasum(pixelEvalArrSize, tempMCSHResult, 1) * dv;
				stencil[index] = sumArr(tempMCSHResult, pixelEvalArrSize) * dv;
				index++;
			}
		}
	}
	free(refX);
	free(refY);
	free(refZ);
	free(tempXArr);
	free(tempYArr);
	free(tempZArr);
	free(tempMCSHResult);*/