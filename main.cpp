#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <thread>
#include <filesystem>
#define STB_IMAGE_IMPLEMENTATION
#include "Dependencies/stb_image.h"
#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))
#define GENE_COUNT 512
#ifndef POPULATION_SIZE
	#define POPULATION_SIZE 15
#endif
//The number of divisions that the sample space should be divided into, this is the power for 2, i.e. value of 2 is 2 ^ 2 divisions in each axis
//E.g. a value of 2 would split the space into 4 square in each axis, making 16 squares in total
//Making this smaller will increase the number of generations performed per seconds but has a trade off in terms of accuracy, i.e. it'll take more generations to get the same result
#ifndef SAMPLE_SPACE_POWER
	#define SAMPLE_SPACE_POWER 4 //256 samples per member
#endif
//The total number of sample being made, this is the number of squares in the sample space
#define SAMPLE_SIZE 1 << (SAMPLE_SPACE_POWER << 1)

static std::chrono::steady_clock::time_point startTime;
static std::string outputPath = "";
static int lastGeneration = -1;	//The generation to end on
//Only used in frames mode
static bool frameMode = false;	//Stores whether the ga is being ran in frames mode
static int frameIndex = 0;	//Stores the index into the frame buffer

//------Structures------
struct Gene {
	//The first 6 bytes of a gene are the positions, the position is treated as -128 to 127 for wraparound
	//The last 4 is the colour
	unsigned char data[10];
	//Stores whether the gene was used in natural selection, if this ends as false, the gene is deleted
	bool geneUsed = false;
};
struct Member {
	//Every 10 bytes is a new gene
	Gene* genes[GENE_COUNT];
	//Stores the fitness of the member
	int fitness = 0;
};
struct Population {
	//Stores a set of members
	Member members[POPULATION_SIZE];// = { 0 };
	//Stores a mating pool
	std::vector<Member*>* matingPool;
	//Stores a reference to the best member
	Member* bestMember;
	//Stores the largest fitness of a member in the population
	int largestFitness = INT_MIN;
	//Stores the smallest fitness of a member in the population
	int smallestFitness = INT_MAX;
};

//------Procedures/Functions------
//------Let there be light
//Creates the genes for a member
//Arguments:	-A member	-Member*
void generateGenes(Member* _member) {
	for (int i = 0; i < GENE_COUNT; i++) {
		_member->genes[i] = new Gene();
		int index = i * 10;
		int newValue = 0;
		_member->genes[i]->data[0] = rand() & 0xFF; //x0
		_member->genes[i]->data[1] = rand() & 0xFF; //y0
		//Ensures that all of the coordinates are within bounds
		newValue = _member->genes[i]->data[0] + rand() % 20; //x1
		_member->genes[i]->data[2] = ((newValue < 0) ? 256 - (newValue & 0xFF) : (newValue >= 256) ? 255 - newValue & 0xFF : newValue);
		newValue = _member->genes[i]->data[1] + rand() % 20; //y1
		_member->genes[i]->data[3] = ((newValue < 0) ? 256 - (newValue & 0xFF) : (newValue >= 256) ? 255 - newValue & 0xFF : newValue);
		newValue = _member->genes[i]->data[0] + rand() % 20; //x2
		_member->genes[i]->data[4] = ((newValue < 0) ? 256 - (newValue & 0xFF) : (newValue >= 256) ? 255 - newValue & 0xFF : newValue);
		newValue = _member->genes[i]->data[1] + rand() % 20; //y3
		_member->genes[i]->data[5] = ((newValue < 0) ? 256 - (newValue & 0xFF) : (newValue >= 256) ? 255 - newValue & 0xFF : newValue);

		//Can use & rather than % as & is check across the whole range of values, not limitng any (if & was using 0b11111110 then the 'random' numbers would only be even)
		_member->genes[i]->data[6] = rand() & 0xFF; //r
		_member->genes[i]->data[7] = rand() & 0xFF; //g
		_member->genes[i]->data[8] = rand() & 0xFF; //b
		_member->genes[i]->data[9] = rand() & 0xFF; //a
	}
}
//Creates a new population
//Arguments:	-A population container	-Population*
void createPopulation(Population* _population) {
	for (int i = 0; i < POPULATION_SIZE; i++) {
		//Generates the genes for the member
		generateGenes(&_population->members[i]);
	}
}
//Mutates a gene
//Arguments:	-The original gene	-Gene*
//Returns:		-A new gene			-Gene*
Gene* mutate(Gene* _gene) {
	//Creates a new gene
	Gene* newGene = new Gene();
	//Mutates all values, including coords and colours
	for (int i = 0; i < 10; i++) {
		int newValue = _gene->data[i];
		//Mutates the value by 3 if it is a coord and by 20 if it is a colour
		newValue += (i < 6) ? (rand() % 7) - 3 : (rand() % 41) - 20;
		//Normalises the data to between 0 and 256
		newGene->data[i] = ((newValue < 0) ? 256 - (newValue & 0xFF) : (newValue >= 256) ? 255 - newValue & 0xFF : newValue);
	}
	//Returns the new gene
	return newGene;
}
//Exports the best member of a population
//Arguments:	-The population		-Population*
//				-The output name	-std::string
//				-The generation		-Int
void exportBestMember(Population* population, std::string name, int generation) {
	std::ofstream wf(outputPath + name + std::string(".dat"), std::ios::out | std::ios::binary);
	char gen[4];//Makes the first 4 bytes the generation number
	gen[0] = (generation >> 24) & 0xFF;
	gen[1] = (generation >> 16) & 0xFF;
	gen[2] = (generation >> 8) & 0xFF;
	gen[3] = (generation >> 0) & 0xFF;
	wf.write(gen, 4);
	for (int i = 0; i < GENE_COUNT; i++) wf.write((char*)population->bestMember->genes[i]->data, 10);
	wf.close();
}

//Loads an image into the targetPixelData buffer
//Arguments:	-The image path			-Char*
//				-The target pixel data	-Unsigned char*
//Returns:		-Success				-Bool
bool loadTargetImage(const char* imagePath, unsigned char* pixelData) {
	//Gets the target image
	int width, height, channels;
	unsigned char* imageData = stbi_load(imagePath, &width, &height, &channels, 3);
	//Requires that the width and height are 256 x 256
	if (width != 256 || height != 256) {
		std::cout << "The image has to be 256 x 256" << std::endl;
		return false;
	}
	//Copies the image data into the buffer
	for (int i = 0; i < 65536 * 3; i++) {
		pixelData[i] = imageData[i];
	}
	free(imageData);
	return true;
}

//Perforns natural selection on the generation
//Arguments:	-The old population -Population
//				-The new population	-Population*
void performNaturalSelection(Population* oldPopulation, Population* newPopulation) {
	//Stores all of the randomly chosen parents, each two members are the parents for a new member
	Member* parents[POPULATION_SIZE*2];
	//Randomly fills the parents array with members, making sure that the first parent is always the best member
	for (int i = 0; i < POPULATION_SIZE; i++) {
		parents[i * 2] = oldPopulation->bestMember;
		parents[i * 2 + 1] = (*oldPopulation->matingPool)[rand() % oldPopulation->matingPool->size()];
	}
	//For each memebr in the population, genes are assigned
	//It assigns the first gene for all members, then the second gene and so on
	for (int i = 0; i < GENE_COUNT; i++) {
		//Loops through every member
		for (int j = 0; j < POPULATION_SIZE; j++) {
			//Randomly chooses the gene from the two parents assigned to this new member
			int memberIndex = (rand() % 2) + j;
			//Randomly mutates the gene
			if (rand() % 100 < 5) {//Mutate
				newPopulation->members[j].genes[i] = mutate(parents[memberIndex]->genes[i]);
			}
			else {
				parents[memberIndex]->genes[i]->geneUsed = true;	//Sets the flag on the gene to say that it got used
				newPopulation->members[j].genes[i] = parents[memberIndex]->genes[i];
			}
		}
		//Stores a buffer of pointers to genes that have been handled
		Gene* handledGenes[POPULATION_SIZE] = { nullptr };
		int handledGenesSize = 0;	//Stores the size of the array
		//Loops through every member of the old population, deleting genes if they weren't used, this is obtained via the geneUsed flag
		for (int j = 0; j < POPULATION_SIZE; j++) {
			//Loops through the buffer to see if the current gene has already been handled
			int iter = 0;
			for (iter = 0; iter < handledGenesSize && oldPopulation->members[j].genes[i] != handledGenes[iter]; iter++);
			//If the iter is bigger than or equal to the size then the gene isn't in the handledGenes list
			if (iter >= handledGenesSize) {
				//Checks if the gene should be deleted
				if (!oldPopulation->members[j].genes[i]->geneUsed) delete oldPopulation->members[j].genes[i];
				else oldPopulation->members[j].genes[i]->geneUsed = false;
				//Adds the gene to handledGenes to signify that a duplicate gene should not be handled
				handledGenes[handledGenesSize] = oldPopulation->members[j].genes[i];
				handledGenesSize++;
			}
		}
	}
}

//------Populous evaluation
//Calculates the fitness of a mamber
//Arguments:	-The population		-Population*
//				-The member			-Member*
//				-The sample space	-Unsigned char*
//				-The target image pixel data	-Unsigned char*
void calculateMemberFitness(Population* _population, Member* _member, unsigned char* sampleSpace, unsigned char* targetImagePixels) {
	//------Calculates the fitness for the member of the generation
	//Stores the colours of the points, every 3 elements is a pixel
	unsigned char render[3 * SAMPLE_SIZE] = { 0 };
	//Loops through every gene, evaluating all the random points for each
	for (int i = 0; i < GENE_COUNT; i++) {
		Gene* gene = _member->genes[i];	//Gets the gene
		unsigned char* verticies = gene->data;
		//Calculates the axis-aligned bounding box of the gene's polygon
		int minx, miny, maxx, maxy;
		minx = min(min(verticies[0], verticies[2]), verticies[4]);
		miny = min(min(verticies[1], verticies[3]), verticies[5]);
		maxx = max(max(verticies[0], verticies[2]), verticies[4]);
		maxy = max(max(verticies[1], verticies[3]), verticies[5]);
		//Calculates the area of the whole gene polygon
		//Uses the fact that the sum of the areas that the point makes with all sides is equal to the area of the whole triangle
		//2A = (x2-x1)(y0-y1) + (y2-y1)(x1-x0)
		int wholeArea = (verticies[4] - verticies[2]) * (verticies[1] - verticies[3]) + (verticies[5] - verticies[3]) * (verticies[2] - verticies[0]);
		//Loops through every random point
		for (int j = 0; j < SAMPLE_SIZE; j++) {
			int x = sampleSpace[j << 1], y = sampleSpace[(j << 1) + 1];	//Gets the x and y of the random point
			//Checks if the point is within the axis aligned bounding box of the gene's polygon
			if (x < minx || x > maxx || y < miny || y > maxy) continue;
			//Calculates the area between the point and the first side
			int area_0 = (x - verticies[2]) * (verticies[1] - verticies[3]) + (y - verticies[3]) * (verticies[2] - verticies[0]);
			//Calculates the area between the point and the second side
			int area_1 = (x - verticies[4]) * (verticies[1] - verticies[5]) + (y - verticies[5]) * (verticies[4] - verticies[0]);
			//Calculates the area between the point and the first side
			int area_2 = (x - verticies[2]) * (verticies[5] - verticies[3]) + (y - verticies[3]) * (verticies[2] - verticies[4]);
			//Checks if the random point is within the gene's polygon
			if (abs(area_0) + abs(area_1) + abs(area_2) != abs(wholeArea)) continue;
			//Draws the polygon colour for that point
			//Performs alpha blending between the colour and the render pixel
			float a = ((float)(gene->data[9]) / 255.0f) * 1.0f; //Multiplies by the blend factor, 1.0f
			float c = 1.0f - a;
			render[j * 3] = (unsigned char)(a * (float)gene->data[6] + c * (float)render[j * 3]);
			render[j * 3 + 1] = (unsigned char)(a * (float)gene->data[7] + c * (float)render[j * 3 + 1]);
			render[j * 3 + 2] = (unsigned char)(a * (float)gene->data[8] + c * (float)render[j * 3 + 2]);
		}
	}
	//Sets the fitness of the member to the total difference between the sample image and the target image
	//This is expected to be normalised later along with the rest of the members
	for (int i = 0; i < SAMPLE_SIZE; i++) {
		//Gets the coordinates, i * 1 and +1 for y
		int x = sampleSpace[i << 1], y = sampleSpace[(i << 1) + 1];
		//Calculates the difference between the member pixel and the target image pixel in all three channels
		//and then adds the magnitude of this value to the fitness
		//Render can be indexed by i as it is only a subset of pixels however targetImagePixels contains all pixels and so x any y has to be converted to linear
		_member->fitness += abs(render[i * 3] - targetImagePixels[((y << 8) + x) * 3]);	//r
		_member->fitness += abs(render[i * 3 + 1] - targetImagePixels[((y << 8) + x) * 3 + 1]);	//g
		_member->fitness += abs(render[i * 3 + 2] - targetImagePixels[((y << 8) + x) * 3 + 2]);	//b
	}
	//Checks if the member has the smallest or biggest fitness
	if (_member->fitness > _population->largestFitness) _population->largestFitness = _member->fitness;
	else if (_member->fitness < _population->smallestFitness) {//This member has a better fitness and so it overwrites bestMember
		_population->bestMember = _member;
		_population->smallestFitness = _member->fitness;
	}
}
//Calculates the fitness of every memeber in the population
//Arguments:	-The population					-Population*
//				-The target image pixel data	-Unsigned char*
void calculateFitnesses(Population* _population, unsigned char* targetImagePixels) {
	//Creates the mating pool for the population
	_population->matingPool = new std::vector<Member*>();

	//------Creates a sample space for all of the members to use
	//Generates random points
	unsigned char* sampleSpace = (unsigned char*)malloc(SAMPLE_SIZE * 2);
	//Old algorithm - 5241 microseconds, new 191 microseconds
	int squareSize = 1 << (8 - SAMPLE_SPACE_POWER);	//The width in pixels of a square in the division
	for (int k = 0; k < SAMPLE_SIZE; k++) { //Loops through every square in the grid
		//Calculates the currently in square in i and j space
		int i = k & ((1 << SAMPLE_SPACE_POWER) - 1);
		int j = k >> SAMPLE_SPACE_POWER;
		//Calculates the coordinates
		int x = squareSize * i + (rand() & (squareSize - 1));
		int y = squareSize * j + (rand() & (squareSize - 1));
		sampleSpace[k << 1] = x;
		sampleSpace[(k << 1) + 1] = y;
	}

	//------Calculates the fitness for each member in the population, giving it the sampleSpace
	//Stores all of the threads for each member
	std::thread* memberThreads[POPULATION_SIZE];
	for (int i = 0; i < POPULATION_SIZE; i++) {
		_population->members[i].fitness = 0;
		std::thread* memberThread = new std::thread(calculateMemberFitness, _population, &_population->members[i], sampleSpace, targetImagePixels);
		memberThreads[i] = memberThread;
		//calculateMemberFitness(_population, &_population->members[i], sampleSpace, targetImagePixels);
	}
	//Joins all the threads
	for (int i = 0; i < POPULATION_SIZE; i++) {
		memberThreads[i]->join();
		delete memberThreads[i];	//Deletes the heap allocated thread
	}

	free(sampleSpace);
	//------Normalises the fitness of each member to between 1 and 10
	//Calculates the gradient to be used in the normalisation equation
	double gradient = 9 / (double)((long)_population->smallestFitness - (long)_population->largestFitness);
	for (int i = 0; i < POPULATION_SIZE; i++) {
		//Calculates the normalised fitness for the member, n = (f - s) * m + 10
		_population->members[i].fitness = (int)((double)(_population->members[i].fitness - _population->smallestFitness) * gradient + 10);

		//Checks if the current member is the best member, it does this by comparing the pointers
		if (&_population->members[i] == _population->bestMember) continue;
		//Adds the member into the mating pool based on its normalised fitness
		for (int j = 0; j < _population->members[i].fitness; j++) {
			_population->matingPool->push_back(&_population->members[i]);
		}
	}
}

//Runs a generation
//Arguments:	-The population					-Population*
//				-The generation number			-Int
//				-The target image pixel data	-unsigned char*
//				-The frame buffer (if any)		-std::vector
//Returns:		-The new population				-Population*
Population* runGeneration(Population* population, int genNumber, unsigned char* targetImagePixels, std::vector<std::pair<int, std::string>>* frameBuffer = nullptr) {
	std::cout << "\r------Generation " << genNumber << ", Time: ";
	double timeElapsed = (double)std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - startTime).count() / (double)1000000;
	std::cout << timeElapsed << "s, Avg: "; 
	std::cout << ((double)genNumber / timeElapsed) << "gens/s------";
	//Calculates the fitness values for each member in the population
	calculateFitnesses(population, targetImagePixels);
	if (!frameMode) {
		//Exports the best member every 50 generations
		if (genNumber % 1000 == 0) exportBestMember(population, "bestMember", genNumber);
		//Exports the best member every power of 2
		if (genNumber != 0 && (genNumber & (genNumber - 1)) == 0) exportBestMember(population, "bestMember_gen" + std::to_string(genNumber), genNumber);
	}
	else {
		//If the ga is in frames mode, it checks to see if the number of generations has been met for the target image to be changed
		if (frameIndex < frameBuffer->size() && (*frameBuffer)[frameIndex].first == genNumber) {
			exportBestMember(population, "frame_" + std::to_string(frameIndex - 1), genNumber);
			std::cout << "\nChanging to frame index: " << (frameIndex - 1) << std::endl;
			//Changes the current frame
			if (!loadTargetImage((*frameBuffer)[frameIndex].second.c_str(), targetImagePixels)) { return 0; }
			frameIndex++;
		}
	}
	//Checks if the generation is the last generaton
	if (lastGeneration == genNumber) exportBestMember(population, (frameMode) ? ("frame_" + std::to_string(frameIndex - 1)) : "lastMember", genNumber);
	//Generates a new population
	Population* newPopulation = new Population();	
	performNaturalSelection(population, newPopulation);
	//Deletes the old population
	delete population->matingPool;
	delete population;

	return newPopulation;
}

int main(int argc, char* argv[]) {
	//Stores the image path
	char* imgPath = (char*)"";
	//Extracts the cli arguments
	for (int i = 0; i < argc; i++) {
		if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) {
			imgPath = argv[i + 1];
			i++;
		}
		else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
			outputPath = std::string(argv[i + 1]);
			i++;
		}
		else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--frames") == 0) {
			frameMode = true;
		}
		else if (strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--last") == 0) {
			lastGeneration = std::stoi(argv[i+1]);
			i++;
		}
		else if (i != 0) {
			std::cout << "Invalid argument: " << argv[i] << std::endl;
		}
	}
	if (*imgPath == 0x00) {
		std::cout << ((!frameMode) ? "Please provide an input image" : "Please provide a frame buffer path") << std::endl;
		return 0;
	}
	if (outputPath.length() == 0) {
		std::cout << "Please provide an output path" << std::endl;
		return 0;
	}
	
	//Stores the pixel data
	unsigned char* targetPixelData = (unsigned char*)malloc(65536 * 3);
	//------If the ga is in frames mode, then it generates an index of the input directory, and uses that as the target images
	//Generates the relationships between generation number and path
	std::vector<std::pair<int, std::string>>* frameBuffer = new std::vector<std::pair<int, std::string>>();
	if (frameMode) {
		for (const auto& entry : std::filesystem::directory_iterator(imgPath)) {
			std::string path = entry.path().string();
			//Stores the index of the start of the number
			int startIndex = 0;
			int endIndex = 0;
			//Reverses through the path list, finding the indicies of the gen number
			for (int i = path.length() - 1; i > 0 && (startIndex == 0 || endIndex == 0); i--) {
				if (path[i] == '.') endIndex = i;
				if (path[i] == '_') startIndex = i;
			}
			if (endIndex - startIndex <= 1) continue;
			char* genNumber = (char*)malloc(endIndex - startIndex);
			memset(genNumber, 0x00, endIndex - startIndex);
			for (int i = 0; i < endIndex - startIndex - 1; i++) {
				genNumber[i] = path[i + startIndex + 1];
			}
			//Adds the frame relationship to the frameBuffer
			std::pair<int, std::string> frameRelationship(std::stoi(genNumber), path);
			frameBuffer->push_back(frameRelationship);

			free(genNumber);
		}
		std::sort(frameBuffer->begin(), frameBuffer->end());
		//Loads the first frame in the frame buffer into the target image
		if (!loadTargetImage((*frameBuffer)[0].second.c_str(), targetPixelData)) { return 0; }
		frameIndex = 1;	
	}
	else {
		//Loads the target image
		if (!loadTargetImage(imgPath, targetPixelData)) { return 0; }
	}

	//Creates the initial population
	Population* population = new Population();
	createPopulation(population);

	//Runs the ga
	startTime = std::chrono::steady_clock::now();
	int genNumber = 0;
	while (true){
		//Runs a generation
		population = runGeneration(population, genNumber, targetPixelData, (frameMode) ? frameBuffer : nullptr);
		if (genNumber == lastGeneration) return 0;	//Exits if the generation is the last generation, defined in the cli
		genNumber++;
	}


	std::string in;
	std::cin >> in;
}

//auto t1 = std::chrono::high_resolution_clock::now();
//auto t2 = std::chrono::high_resolution_clock::now();
//auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
//std::cout << duration << std::endl;


//std::ofstream wf("C:/Users/Aidan/Desktop/out.dat", std::ios::out | std::ios::binary);
//wf.write((char*)randPoints, SAMPLE_SIZE * 2);
//wf.write((char*)render, SAMPLE_SIZE * 3);
//wf.close();
