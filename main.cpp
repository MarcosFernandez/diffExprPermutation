/*
 * main.cpp
 *
 *  Created on: 30 Xan, 2015
 *      Author: marcos
 */
#include "common.h"

/**
 * \brief Print application help
 */
void printHelpPermutations()
{
	cout << "diffExprPermutation: Find differentially expressed genes from a set of control and cases samples using a permutation strategy." << endl;
	cout << "diffExprPermutation -m PERM -f input [--c1 condition1 --c2 condition2 -n permutations -H stopHits -s statistic] -o outFile" << endl;

	cout << "Inputs:" << endl;
	cout << " -m PERM         PERMUTATIONS MODE run the permutations over a set of genes and samples." << endl;
	cout << " -f input        Space separated table. Format: sampleName group lib.size norm.factors gene1 gene2 ... geneN" << endl;
	cout << "Outputs:" << endl;
	cout << " -o outFile      Output file name" << endl;
	cout << "Options:" << endl;
	cout << " --c1	condition1            Condition that determine one of two groups [default: case]" << endl;
	cout << " --c2	condition2            Condition that determine other group [default: control]" << endl;
	cout << " -n PERMUTATIONS             Maximum number of permutations [default: 1000000]"<< endl;
	cout << " -H STOP_HITS                Maximum hits to stop permutations [default: 5]" << endl;
	cout << " -s statistic                Statistic to compute pvalue median|perc25|perc75 [default: median]" << endl;
}


/**
 * \brief Print application help
 */
void printHelpSplitter()
{
	cout << "diffExprPermutation: Splits a file of gene expressions values to perform a parallel running." << endl;
	cout << "diffExprPermutation -m SPLIT -f input -o outPrefix [--n chunks --c1 condition1 --c2 condition2]" << endl;

	cout << "Inputs:" << endl;
	cout << " -m SPLIT        SPLITTER MODE Splits original input file in a set of chunks. The goal is to perform a parallel running of the permutations process." << endl;
	cout << " -f input        Space separated table. Format: sampleName group lib.size norm.factors gene1 gene2 ... geneN" << endl;
	cout << "Outputs:" << endl;
	cout << " -o outPrefix    Prefix for the output mode" << endl;
	cout << "Options:" << endl;
	cout << " --n	chunks            Number of chunks to divide the file [default: 500]" << endl;
	cout << " --c1	condition1        Condition that determine one of two groups [default: case]" << endl;
	cout << " --c2	condition2        Condition that determine other group [default: control]" << endl;
}

/**
 * \brief Print application help
 */
void printHelpFDR()
{
	cout << "diffExprPermutation: Adjust pvalues by False Discovery Rate method from a set of genes distributed in files." << endl;
	cout << "diffExprPermutation -m FDR -d directory -o outFile " << endl;

	cout << "Inputs:" << endl;
	cout << " -m FDR        FDR MODE. P values adjustment per Benjamini-Hochber." << endl;
	cout << " -d directory Directory where are located the files to adjust the pvalue by Benjamini-Hochber method." << endl;
	cout << "               File Format: gene [tab] diff_median [tab] medianCase [tab] medianControl [tab] fold_change [tab] hits [tab] N_Perm [tab] median_pv [tab] median_pv_fdr" << endl;
	cout << "Outputs:" << endl;
	cout << " -o outFile    Output file with FDR corrected values" << endl;
}


/**
 * \brief Check Splitter Arguments
 * \param string fileInput - Name input file
 * \param string outPrefix -  Prefix for output files
 * \param unsigned int chunks - Number of chunks to create
 * \param string condition1 - First condition group. Usually case.
 * \param string condition1 - Second condition group. Usually control.
 */
inline bool argumentCheckingSplitter(const string & fileInput,const string & outPrefix, unsigned int chunks, const string & condition1,const string & condition2)
{
	bool bWrong = false;

	if (fileInput.empty())
	{
		cout << "Sorry!! No input file was specified!!" << endl;
		return true;
	}

	if (outPrefix.empty())
	{
		cout << "Sorry!! No output prefix was specified!!" << endl;
		return true;
	}

	if (chunks > 2000)
	{
		cout << "Sorry!! Too match chunks to create, please select a number under 2000!!" << endl;
		return true;
	}

	if (chunks < 1)
	{
		cout << "Sorry!! Too many chunks to create, please selecte a unmber between 1 and 2000!!" << endl;
		return true;
	}

	if (condition1.empty())
	{
		cout << "Sorry!! Condition group 1 is empty!!" << endl;
		return true;
	}

	if (condition2.empty())
	{
		cout << "Sorry!! Condition group 2 is empty!!" << endl;
		return true;
	}

	return bWrong;
}

/**
 * \brief Check Permutation Arguments
 * \param string fileInput - Name input file
 * \param string outFile -   Name output file
 * \param string condition1 - First condition group. Usually case.
 * \param string condition1 - Second condition group. Usually control.
 * \param string statistics - Statistic selected
 */
inline bool argumentCheckingPermutations(const string & fileInput,const string & outFile, const string & condition1,const string & condition2,const string & statistic)
{
	bool bWrong = false;

	if (fileInput.empty())
	{
		cout << "Sorry!! No input file was specified!!" << endl;
		return true;
	}

	if (outFile.empty())
	{
		cout << "Sorry!! No output file was specified!!" << endl;
		return true;
	}

	if (condition1.empty())
	{
		cout << "Sorry!! Condition group 1 is empty!!" << endl;
		return true;
	}

	if (condition2.empty())
	{
		cout << "Sorry!! Condition group 2 is empty!!" << endl;
		return true;
	}

	if (statistic.compare("median") != 0 && statistic.compare("perc25") && statistic.compare("perc75"))
	{
			cout << "Sorry!! Statistical " << statistic << "is not a valid one!!" << endl;
			return true;
	}

	return bWrong;
}

/**
 * \brief Check FDR Arguments
 * \param string directory - Name input directory
 * \param string outFile -   Name output file
 */
bool argumentCheckingFDR(const string & directory,const string & outFile)
{
	if (directory.empty())
	{
		cout << "Sorry!! No input directory was specified!!" << endl;
		return true;
	}

	if (outFile.empty())
	{
		cout << "Sorry!! No output file was specified!!" << endl;
		return true;
	}

	return false;
}

/**
 * \brief Run Splitter Creates chunks of genes to be run in parallel
 * \param argc Number of arguments
 * \param argv List of arguments
 */
void runSplitter(int argc, char ** argv)
{
	string fileInput;
	string outputPrefix;
	unsigned int chunks = 500;
	string condition1 = "case";
	string condition2 = "control";

	//1. Process Splitter Parameters
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
			fileInput = argv[i+1];
		}
		else if (strcmp(argv[i],"-o") == 0)
		{
			outputPrefix = argv[i+1];
		}
		else if (strcmp(argv[i],"--n") == 0)
		{
			chunks = atoi(argv[i+1]);
		}
		else if (strcmp(argv[i],"--c1") == 0)
		{
			condition1 = argv[i+1];
		}
		else if (strcmp(argv[i],"--c2") == 0)
		{
			condition2 = argv[i+1];
		}
	}

	//4.Checking arguments
	if(argumentCheckingSplitter(fileInput,outputPrefix,chunks,condition1,condition2))
	{
		return;
	}

	//3. Get file header to count number of lines
	string header = getFileHeader(fileInput);
	vector <string> headerFields;
	split(header,headerFields);
	unsigned int numberGenes = headerFields.size() - 4;
	unsigned int genesPerFile = (int) round((float)numberGenes / (float) chunks);

	//4. Load file in memory
	vector <Gene> vGenes; //vector of genes where each gene has a vector of sampleGenes, each sampleGene contains sample name, expression value and group
	//Parsing Input file
	if (!loadFileInfo(fileInput,vGenes,condition1,condition2))
	{
		cout << "Sorry!! Can not open file " << fileInput << endl;
		return;
	}


	//5.Print chunk
	unsigned int genesPrintedInChunk = 0;
	unsigned int genesPrinted = 0;
	unsigned int chuncksCreated = 0;

	vector <string> chunckContents; // vector of lines to be printed in the chunked file

	for (vector<Gene>::const_iterator iter = vGenes.begin(); iter != vGenes.end(); ++iter)
	{
		if (genesPrintedInChunk == 0)
		{
			//Add header contents
			chunckContents.push_back(headerFields[0] + " " + headerFields[1] + " " + headerFields[2] + " " + headerFields[3]);
			//Add Gene
			addGeneContents(true,chunckContents, (*iter), condition1, condition2);
			//Extreme case just one gene per file
			if(genesPerFile == 1)
			{
				genesPrinted ++;
				chuncksCreated ++;
				//print list of genes to file chunk
				printChunck(chunckContents,chuncksCreated,outputPrefix);
				//Reset counters
				genesPrintedInChunk = 0;
			}
			else
			{
				genesPrintedInChunk ++;
			}
		}
		else if (genesPrintedInChunk == genesPerFile -1)
		{
			//Last gene to add
			addGeneContents(false,chunckContents, (*iter), condition1, condition2);
			genesPrinted += genesPrintedInChunk + 1;
			chuncksCreated ++;
			//print list of genes to file chunk
			printChunck(chunckContents,chuncksCreated,outputPrefix);
			//Reset counters
			genesPrintedInChunk = 0;
		}
		else
		{
			//Add new Gene
			addGeneContents(false,chunckContents, (*iter), condition1, condition2);
			genesPrintedInChunk ++;
		}
	}

	//Remaining not printed
	if(!chunckContents.empty())
	{
		genesPrinted += genesPrintedInChunk;
		chuncksCreated ++;
		//print list of genes to file chunk
		printChunck(chunckContents,chuncksCreated,outputPrefix);
	}

	cout << "Split finished!!" << endl;
	cout << " Genes per file: " << genesPerFile << endl;
	cout << " Genes printed: " << genesPrinted << endl;
	cout << " Chunks Created: " << chuncksCreated << endl;
}

/**
 * \brief Run Permutations Run the applications performing permutation to calculate a p value
 * \param argc Number of arguments
 * \param argv List of arguments
 */
void runPermutations(int argc, char ** argv)
{
	string fileInput = "";
	string outFile = "";
	string condition1 = "case";
	string condition2 = "control";
	unsigned int nPermutations = 1000000;
	unsigned int stopHits = 5;
	string statistic = "median";
	double fStatisticValue = 0;
	bool doMedian = true;

	vector <Gene> vGenes; //vector of genes where each gene has a vector of sampleGenes, each sampleGene contains sample name, expression value and group
	/**
	 * BRACA1 -> A,true,0.75
	 *        -> B,false,0.85
	 *        ...
	 * BRACA2 -> A,true,0.15
	 *        -> B,false,0.20
	 *        ...
	 */

	//1.Process PERMUTATION parameters
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
			fileInput = argv[i+1];
		}
		else if (strcmp(argv[i],"-o") == 0)
		{
			outFile = argv[i+1];
		}
		else if (strcmp(argv[i],"--c1") == 0)
		{
			condition1 = argv[i+1];
		}
		else if (strcmp(argv[i],"--c2") == 0)
		{
			condition2 = argv[i+1];
		}
		else if (strcmp(argv[i],"-n") == 0)
		{
			nPermutations = atoi(argv[i+1]);
		}
		else if (strcmp(argv[i],"-H") == 0)
		{
			stopHits = atoi(argv[i+1]);
		}
		else if (strcmp(argv[i],"-s") == 0)
		{
			statistic = argv[i+1];
		}
	}

	//Check Arguments
	if(argumentCheckingPermutations(fileInput,outFile,condition1,condition2,statistic))
	{
		return;
	}

	//Updates statistic
	string headerOutput = "gene\tdiff_median\tmedianCase\tmedianControl\tfold_change\thits\tN_Perm\tmedian_pv\tmedian_pv_fdr";

	if	(statistic.compare("perc25") == 0)
	{
		fStatisticValue = 0.25;
		doMedian = false;
		headerOutput = "gene\tdiff_LowQ\tmedianCase\tmedianControl\tfold_change\thits\tN_Perm\tlowerq_pv\tlowerq_pv_fdr";
	}
	else if	(statistic.compare("perc75") == 0)
	{
		fStatisticValue = 0.75;
		doMedian = false;
		headerOutput = "gene\tdiff_UpQ\tmedianCase\tmedianControl\tfold_change\thits\tN_Perm\tupperq_pv\tupper_pv_fdr";
	}

	//Parsing Input file
	if (!loadFileInfo(fileInput,vGenes,condition1,condition2))
	{
		cout << "Sorry!! Can not open file " << fileInput << endl;
		return;
	}

	//Montecarlo permutations per each Gene
	vector <double> vPvalues;

	for (vector<Gene>::iterator iter = vGenes.begin(); iter != vGenes.end(); ++iter)
	{

		//1. Get Difference median control and cases
		if(doMedian)
		{
			(*iter).originalDiff = (*iter).diffMedian();
		}
		else
		{
			(*iter).originalDiff = (*iter).diffPercentile(fStatisticValue);
		}

		//Permutation loop
		double permutationDiference = 0;
		for (unsigned int i = 0; i < nPermutations; i++)
		{
			//Try to stop if the limit of hits is reached
			if ((*iter).hits == stopHits)
			{
				(*iter).pValue = ((double)(*iter).hits / (double) (i));
				(*iter).nPermutationsDone = i;
				break;
			}
			//Random shuffle of values
			(*iter).random_shuffle();
			//Calculates new difference for the new set of groups
			if(doMedian)
			{
				permutationDiference = (*iter).diffMedian();
			}
			else
			{
				permutationDiference = (*iter).diffPercentile(fStatisticValue);
			}

			//Checks for hits
			if(permutationDiference > (*iter).originalDiff)
			{
				(*iter).hits ++;
			}
		}

		//Those genes that its number of hits is below stopHits value estimate its pValue
		if ((*iter).pValue < 0)
		{
			(*iter).pValue = ((double)(*iter).hits + 1) / (double) nPermutations;
			(*iter).nPermutationsDone = nPermutations;
		}

		//Add pvalue to a vector of pvalues for being correcte by FDR
		vPvalues.push_back((*iter).pValue);
	}

	vector <double> correctedPvalues;
	correct_pvalues_fdr(vPvalues,correctedPvalues);

	//Print to file
	std::ofstream outfile (outFile.c_str(),std::ofstream::out);


	//Header File
	outfile << headerOutput << endl;

	vector<double>::const_iterator iterCorrected = correctedPvalues.begin();
	outfile.precision(15);

	for (vector<Gene>::const_iterator iter = vGenes.begin(); iter != vGenes.end(); ++iter)
	{

		outfile << (*iter).geneName << "\t" << (*iter).originalDiff << "\t" << (*iter).originalMedianCases;
		outfile << "\t" << (*iter).originalMedianControl << "\t" << (*iter).foldChange << "\t" << (*iter).hits;
		outfile << "\t" << (*iter).nPermutationsDone << "\t" << (*iter).pValue << "\t" << (*iterCorrected) << endl;

		++ iterCorrected;
	}

	outfile.close();
}


/**
 * \brief Calculates the false discovery rate from a set of gene pvalues
 * \param argc Number of arguments
 * \param argv List of arguments
 */
void runFDR(int argc, char ** argv)
{
	string directory;
	string outFile;

	//1.Process PERMUTATION parameters
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i],"-d") == 0)
		{
			directory = argv[i+1];
		}
		else if (strcmp(argv[i],"-o") == 0)
		{
			outFile = argv[i+1];
		}
	}

	//2. Check Arguments
	if(argumentCheckingFDR(directory,outFile))
	{
		return;
	}

	//3. Read files from folder
	vector <GenePvalue> genePvalues;
	string header;
	loadGenePValuesFromFolder(genePvalues,directory,header);

	//4. Get all pvalues
	vector <double> vPvalues;

	for (vector<GenePvalue>::iterator iter = genePvalues.begin(); iter != genePvalues.end(); ++iter)
	{
		vPvalues.push_back((*iter).getPv());
	}

	//5. Calculate adjusted pvalues
	vector <double> correctedPvalues;
	correct_pvalues_fdr(vPvalues,correctedPvalues);

	//6. Update pvalues
	vector<double>::iterator pvAdjusted = correctedPvalues.begin();
	for (vector<GenePvalue>::iterator iter = genePvalues.begin(); iter != genePvalues.end(); ++iter)
	{
		(*iter).setPvAdjusted((*pvAdjusted));
		++ pvAdjusted;
	}

	//7. Print to output file
	printCorrectedPvalues(genePvalues,header, outFile);
}

int main (int argc, char ** argv)
{
	unsigned int runMode = NO_MODE;

	//1 Parse Arguments to get MODE
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i],"-m") == 0)
		{
			if( strcmp(argv[i+1],"SPLIT") == 0)
			{
				runMode = MODE_SPLITTER;
			}
			else if( strcmp(argv[i+1],"PERM") == 0)
			{
				runMode = MODE_PERMUTATIONS;
			}
			else if( strcmp(argv[i+1],"FDR") == 0)
			{
				runMode = MODE_FDR;
			}
		}
		else if (strcmp(argv[i],"-h") == 0)
		{
			switch (runMode)
			{
				case MODE_SPLITTER:
					printHelpSplitter();
					break;
				case MODE_PERMUTATIONS:
					printHelpPermutations();
					break;
				case MODE_FDR:
					printHelpFDR();
					break;
				default:
					printHelp();
					break;
			}
			return 0;
		}
	}

	//Run the appropriate mode
	switch (runMode)
	{
		case MODE_SPLITTER:
			runSplitter(argc,argv);
			break;
		case MODE_PERMUTATIONS:
			runPermutations(argc,argv);
			break;
		case MODE_FDR:
			runFDR(argc,argv);
			break;
		default:
			printHelp();
			break;
	}


	return 0;
}
