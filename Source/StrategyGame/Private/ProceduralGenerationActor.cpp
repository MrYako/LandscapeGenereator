// Fill out your copyright notice in the Description page of Project Settings.


#include "ProceduralGenerationActor.h"

// Sets default values
AProceduralGenerationActor::AProceduralGenerationActor()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = false;
	size = 64;
	roughness = 1.0;
}

// Called when the game starts or when spawned
void AProceduralGenerationActor::BeginPlay()
{
	Super::BeginPlay();
	
}

// Called every frame
void AProceduralGenerationActor::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

TArray<FVector>  AProceduralGenerationActor::DiamondSquareAlgh(int sizeMatrix, float roughnessCoef, bool powHeigh)
{
	heighmap.Empty();
	coordVert.Empty();
	roughness = roughnessCoef;
	size = sizeMatrix + 1; // DiamondSq have 2^n + 1 vertix
	int i, j;
	heighmap.Init(0, size*size);

	heighmap[0] = FMath::RandRange(0.3f, 0.6f);
	heighmap[size * (size - 1)] = FMath::RandRange(0.3f, 0.6f);
	heighmap[size * (size - 1) + size - 1] = FMath::RandRange(0.3f, 0.6f);
	heighmap[size - 1] = FMath::RandRange(0.3f, 0.6f);
	//UE_LOG(LogTemp, Warning, TEXT("%f %f %f %f"), heighmap[0], heighmap[size * (size - 1)], heighmap[size * (size - 1) + size - 1], heighmap[size - 1]);
	
	for (int l = (size - 1); l > 0; l /= 2)
		for (int x = 0; x < size - 1; x += l)
			for (int y = 0; y < size - 1; y += l)
				DiamondSquare(x, y, x + l, y + l);

	FVector currentPoint;
	for (i = 0; i < sizeMatrix; i++)
		for (j = 0; j < sizeMatrix; j++) {
			currentPoint.X = i;
			currentPoint.X = currentPoint.X / sizeMatrix;
			currentPoint.Y = j;
			currentPoint.Y = currentPoint.Y / sizeMatrix;
			if (powHeigh)
				currentPoint.Z = pow( heighmap[i * size + j],2);
			else
				currentPoint.Z = abs(heighmap[i * size + j]);
			coordVert.Add(currentPoint);
		}
	return(coordVert);
}


void  AProceduralGenerationActor::DiamondSquare(int lx, int ly, int rx, int ry)
{
	int l = (rx - lx) / 2;
	
	
	if (l > 0)
	{
		//UE_LOG(LogTemp, Warning, TEXT("l = %d lx = %d ly = %d rx = %d ry = %d"), l, lx,ly,rx,ry);
		Square(lx, ly, rx, ry);

		Diamond(lx, ly + l, l);
		Diamond(rx, ry - l, l);
		Diamond(rx - l, ry, l);
		Diamond(lx + l, ly, l);
	}
}

void AProceduralGenerationActor::Square(int lx, int ly, int rx, int ry)
{
	int l = (rx - lx) / 2;
	//UE_LOG(LogTemp, Warning, TEXT(" ll =  %d"), l);
	float a = heighmap[size * lx + ly];              //  B--------C
	float b = heighmap[size * lx + ry];              //  |        |
	float c = heighmap[size * rx + ry];              //  |   ce   |
	float d = heighmap[size * rx + ly];              //  |        |
											         //  A--------D
	int cex = lx +l;
	int cey = ly +l;

	heighmap[size * cex + cey] = (a + b + c + d) / 4 + FMath::RandRange(-l * 2 * roughness / size, l * 2 * roughness / size);
	//UE_LOG(LogTemp, Warning, TEXT("heighmap[size * cex + cey] = %f %d"), heighmap[size * cex + cey], size * cex + cey);
}

void AProceduralGenerationActor::Diamond(int tgx, int tgy, int l)
{
	float a, b, c, d;

	if (tgy - l >= 0)
		a = heighmap[size * tgx + tgy - l];
	else
		a = 0;// heighmap[size * tgx + size - l];

	if (tgx - l >= 0)
		b = heighmap[size * (tgx - l) + tgy];
	else
		b = 0; // heighmap[size * (size - l) + tgy];

	if (tgy + l <= size - 1)
		c = heighmap[size * tgx + tgy + l];
	else
		c = 0; // heighmap[size * tgx + l];

	if (tgx + l <= size - 1)
		d = heighmap[size * (tgx + l) + tgy];
	else
		d = 0;  //heighmap[size * l + tgy];

	heighmap[size * tgx + tgy] = (a + b + c + d) / 4 + FMath::RandRange(-l * 2 * roughness / size, l * 2 * roughness / size);
}

TArray<FVector>  AProceduralGenerationActor::CombinationAlgh(float coef1, const TArray<FVector>& coordVector1,float coef2, const TArray<FVector>& coordVector2, int sizeMatrix)
{
	TArray<FVector> result;
	FVector currentPoint;
	int i, j;
	for (i=0;i<sizeMatrix; i++)
		for (j = 0; j < sizeMatrix; j++) {
			currentPoint.X = coordVector1[i * sizeMatrix + j].X;
			currentPoint.Y = coordVector1[i * sizeMatrix + j].Y;
			currentPoint.Z = coordVector1[i * sizeMatrix + j].Z * coef1 + coordVector2[i * sizeMatrix + j].Z * coef2;
			result.Add(currentPoint);
		}
	return result;
}

void AProceduralGenerationActor::thermalErosion(float talus, TArray<FVector> ArrayVert, TArray<FVector> &ArrayVertOut) {
	float coef = 0.2f;
	int sizeMatrix = sqrt(ArrayVert.Num());
	int i, j, k;
	int cur;
	float dtotal;
	float dmax;
	FVector2D di;
	TArray<FVector2D> diArray;
	TArray<int> neighbourhood;

	for (i = 0; i < sizeMatrix; i++)
	{
		for (j = 0; j < sizeMatrix; j++)
		{
			diArray.Empty();
			cur = i * sizeMatrix + j;// current cell
			dtotal = 0;
			dmax = 0;

			//neighbourhood = VonNeimanNeighbourhood(i, j, sizeMatrix);
			//neighbourhood = ModVonNeimanNeighbourhood(i, j, sizeMatrix);
			neighbourhood = MoreNeighbourhood(i, j, sizeMatrix);
			
			for (k = 0; k < neighbourhood.Num(); k++) {
				di.X = neighbourhood[k];
				di.Y = ArrayVert[cur].Z - ArrayVert[di.X].Z;
				if (di.Y > talus) {
					diArray.Add(di);
					dtotal += di.Y;
					if (di.Y > dmax)
						dmax = di.Y;
				}
			}

			if (dmax) {
				ArrayVert[cur].Z = ArrayVert[cur].Z - coef * (dmax - talus);
			}
			for (k = 0; k < diArray.Num(); k++) {
				ArrayVert[diArray[k].X].Z += coef * (dmax - talus) * diArray[k].Y / dtotal;
			}
		}
	}
	ArrayVertOut = ArrayVert;
}

void AProceduralGenerationActor::onlyMinThermalErosion(float talus, TArray<FVector> ArrayVert, TArray<FVector>& ArrayVertOut) {
	float coef = 0.5f;
	int sizeMatrix = sqrt(ArrayVert.Num());
	int i, j, k;
	int cur;
	FVector2D dmax;
	FVector2D di;
	TArray<int> neighbourhood;

	for (i = 0; i < sizeMatrix; i++)
	{
		for (j = 0; j < sizeMatrix; j++)
		{
			cur = i * sizeMatrix + j;// current cell
			dmax.X = -1;
			dmax.Y = talus;

			//neighbourhood = VonNeimanNeighbourhood(i, j, sizeMatrix);
			neighbourhood = ModVonNeimanNeighbourhood(i, j, sizeMatrix);
			//neighbourhood = MoreNeighbourhood(i, j, sizeMatrix);

			for (k = 0; k < neighbourhood.Num(); k++) {
				di.X = neighbourhood[k];
				di.Y = ArrayVert[cur].Z - ArrayVert[di.X].Z;
				if (di.Y > dmax.Y) {
					dmax.X = di.X;
					dmax.Y = di.Y;
				}
			}

			if (dmax.X !=-1 ) {
				ArrayVert[cur].Z = ArrayVert[cur].Z - coef * (dmax.Y - talus);
				ArrayVert[dmax.X].Z += coef * (dmax.Y - talus);
			}
		}
	}
	ArrayVertOut = ArrayVert;
}



void AProceduralGenerationActor::hydraulicErosion(TArray<FVector> ArrayVert, TArray<FVector>& ArrayVertOut, int CountIter, TArray<FString>& Data, bool OptimizationVerson,  bool WriteData) {
	TArray<float> water;
	TArray<float> material;
	material.AddZeroed(ArrayVert.Num());
	water.AddZeroed(ArrayVert.Num());
	int sizeMatrix = sqrt(ArrayVert.Num());
	int i, j;
	float Kr = 0.05f; // count water every iteration 
	float Ks = 0.05f; // is the solubility constant of the terrain
	float Ke = 0.2f; // evaporation coefficient
	float Kc = 0.01f; //  sediment capacity coefficient
	Data.Empty();


	for (j = 0; j < CountIter; j++) {

		//step 1 - add water
		for (i = 0; i < ArrayVert.Num(); i++) {
			water[i] += Kr;
		}

		//step 2 - water eroding 
		for (i = 0; i < ArrayVert.Num(); i++) {

			if (ArrayVert[i].Z - Ks * water[i] < 0) {
				material[i] = material[i] + ArrayVert[i].Z;
				ArrayVert[i].Z = 0;
			}
			else {
				ArrayVert[i].Z = ArrayVert[i].Z - Ks * water[i];
				material[i] = material[i] + Ks * water[i];
			}
			
		}


		//step 3 - transport water
		if (OptimizationVerson)
			onlyMinTransportWater(material, water, ArrayVert);
		else
			transportWater(material, water, ArrayVert);
		

		//step 4 - Evaporation of water
		float mmax;
		float deltaMaterial;
		for (i = 0; i < ArrayVert.Num(); i++) {
			water[i] = water[i] * (1 - Ke);
			mmax = Kc * water[i];
			deltaMaterial = (0 > material[i] - mmax) ? 0 : material[i] - mmax;
			material[i] = material[i] - deltaMaterial;
			ArrayVert[i].Z += deltaMaterial;
			//UE_LOG(LogTemp, Warning, TEXT("%f id = %d"), deltaMaterial, i);
		}
		if (WriteData) {
			FString curStr;
			curStr = FString::FromInt(j) + ";" + FString::SanitizeFloat(erosionScore(ArrayVert)) + "\n";
			Data.Add(curStr);
		}
	}


	ArrayVertOut = ArrayVert;
}

float   AProceduralGenerationActor::erosionScore(const TArray<FVector>& ArrayVert)
{
	TArray<float> matrixS;
	matrixS.Empty();
	int sizeMatrix = sqrt(ArrayVert.Num());
	int i, j;
	int cur = 0;
	float maxDif; // max differance
	float summ = 0;
	for (i = 0; i < sizeMatrix; i++)
	{
		for (j = 0; j < sizeMatrix; j++)
		{
			cur = i * sizeMatrix + j;// current cell
			//find max differance
			maxDif = 0;
			if ((i - 1 >= 0) && (maxDif < abs(ArrayVert[cur].Z - ArrayVert[cur - sizeMatrix].Z)))// h_(i-1,j)
				maxDif = abs(ArrayVert[cur].Z - ArrayVert[cur - sizeMatrix].Z);
			if ((i + 1 < sizeMatrix) && (maxDif < abs(ArrayVert[cur].Z - ArrayVert[cur + sizeMatrix].Z)))// h_(i+1,j)
				maxDif = abs(ArrayVert[cur].Z - ArrayVert[cur + sizeMatrix].Z);
			if ((j - 1 >= 0) && (maxDif < abs(ArrayVert[cur].Z - ArrayVert[cur - 1].Z)))// h_(i,j-1)
				maxDif = abs(ArrayVert[cur].Z - ArrayVert[cur - 1].Z);
			if ((j + 1 < sizeMatrix) && (maxDif < abs(ArrayVert[cur].Z - ArrayVert[cur + 1].Z)))// h_(i,j+1)
				maxDif = abs(ArrayVert[cur].Z - ArrayVert[cur + 1].Z);
			matrixS.Add(maxDif);
			summ += maxDif;

		}
	}
	float amountS = summ / pow(sizeMatrix, 2);
	float dS = 0;
	for (i = 0; i < matrixS.Num(); i++)
	{
		dS += pow(matrixS[i] - amountS, 2);
	}
	dS = sqrt(dS / pow(sizeMatrix, 2));

	float eScore = dS / amountS;

	return eScore;
}

void AProceduralGenerationActor::transportWater(TArray<float> &material, TArray<float> &water, const TArray<FVector>& ArrayVert) {
	float aAverage;
	float dtotal;
	FVector2D di;
	TArray<FVector2D> diArray;
	int cur;
	float coefWater;
	TArray<float> deltaWater;
	TArray<float> deltaMaterialArr;
	TArray<int> neighbourhood;
	int sizeMatrix = sqrt(ArrayVert.Num());
	int i, j, k;

	for (i = 0; i < sizeMatrix; i++) {
		for (j = 0; j < sizeMatrix; j++) {
			diArray.Empty();
			deltaWater.Empty();
			deltaMaterialArr.Empty();
			cur = i * sizeMatrix + j;// current cell
			aAverage = ArrayVert[cur].Z + water[cur];
			dtotal = 0;
			coefWater = 0;

			neighbourhood = MoreNeighbourhood(i, j, sizeMatrix);

			for (k = 0; k < neighbourhood.Num(); k++) {
				di.X = neighbourhood[k];
				di.Y = ArrayVert[cur].Z + water[cur] - (ArrayVert[di.X].Z + water[di.X]);
				if (di.Y > 0) {
					diArray.Add(di);
					dtotal += di.Y;
				}
			}
			if (dtotal) {
				for (k = 0; k < diArray.Num(); k++) {
					aAverage += ArrayVert[diArray[k].X].Z + water[diArray[k].X];
				}
				aAverage = aAverage / diArray.Num();

				if (aAverage > water[cur])
					coefWater = water[cur];
				else
					coefWater = aAverage;
			}

			for (k = 0; k < diArray.Num(); k++) {
				deltaWater.Add(coefWater * diArray[k].Y / dtotal);
				deltaMaterialArr.Add(material[cur] * deltaWater.Last() / water[cur]);
				water[diArray[k].X] += deltaWater.Last();
				material[diArray[k].X] += deltaMaterialArr.Last();
			}
			material[cur] = material[cur] - material[cur] * coefWater / water[cur];
			water[cur] = water[cur] - coefWater;
		}
	}
}

void AProceduralGenerationActor::onlyMinTransportWater(TArray<float>& material, TArray<float>& water, const TArray<FVector>& ArrayVert) {
	FVector2D di;
	FVector2D dmax;
	TArray<int> neighbourhood;
	int cur;
	float coefWater;
	int sizeMatrix = sqrt(ArrayVert.Num());
	int i, j, k;

	for (i = 0; i < sizeMatrix; i++) {
		for (j = 0; j < sizeMatrix; j++) {
			
			dmax.X = -1;
			dmax.Y = 0;

			cur = i * sizeMatrix + j;// current cell
			coefWater = 0;

			//neighbourhood = MoreNeighbourhood(i, j, sizeMatrix);
			neighbourhood = ModVonNeimanNeighbourhood(i, j, sizeMatrix);

			for (k = 0; k < neighbourhood.Num(); k++) {
				di.X = neighbourhood[k];
				di.Y = ArrayVert[cur].Z + water[cur] - (ArrayVert[di.X].Z + water[di.X]);
				if (di.Y > dmax.Y) {
					dmax.X = di.X;
					dmax.Y = di.Y;
				}
			}
			if (dmax.Y) {
				if ((ArrayVert[dmax.X].Z + water[dmax.X] + ArrayVert[cur].Z + water[cur])/2 - (ArrayVert[dmax.X].Z + water[dmax.X]) > water[cur])
					coefWater = water[cur];
				else
					coefWater = (ArrayVert[dmax.X].Z + water[dmax.X] + ArrayVert[cur].Z + water[cur]) / 2 - (ArrayVert[dmax.X].Z + water[dmax.X]);
			
				water[dmax.X] += coefWater;
				material[dmax.X] += material[cur] * coefWater / water[cur];
			}
			material[cur] = material[cur] - material[cur] * coefWater / water[cur];
			water[cur] = water[cur] - coefWater;
		}
	}
}

TArray<int> AProceduralGenerationActor::VonNeimanNeighbourhood(int i, int j, int sizeMatrix) {
	TArray<int> neighbourhood;
	int cur = i * sizeMatrix + j;
	if (i - 1 >= 0)
		neighbourhood.Add(cur - sizeMatrix);
	if (i + 1 < sizeMatrix) 
		neighbourhood.Add(cur + sizeMatrix);
	if (j - 1 >= 0) 
		neighbourhood.Add(cur - 1);
	if (j + 1 < sizeMatrix) 
		neighbourhood.Add(cur + 1);
	return neighbourhood;
}

TArray<int> AProceduralGenerationActor::ModVonNeimanNeighbourhood(int i, int j, int sizeMatrix) {
	TArray<int> neighbourhood;
	int cur = i * sizeMatrix + j;
	if ((i - 1 >= 0) && ( j - 1 >= 0))
		neighbourhood.Add(cur - sizeMatrix - 1);
	if ((i - 1 >= 0) && (j + 1 < sizeMatrix))
		neighbourhood.Add(cur - sizeMatrix + 1);
	if ((i + 1 < sizeMatrix) && (j - 1 >= 0))
		neighbourhood.Add(cur + sizeMatrix - 1);
	if ((i + 1 < sizeMatrix) && (j + 1 < sizeMatrix))
		neighbourhood.Add(cur + sizeMatrix + 1);
	return neighbourhood;
}

TArray<int> AProceduralGenerationActor::MoreNeighbourhood(int i, int j, int sizeMatrix) {
	TArray<int> neighbourhood;
	neighbourhood.Empty();
	neighbourhood.Append(VonNeimanNeighbourhood(i, j, sizeMatrix));
	neighbourhood.Append(ModVonNeimanNeighbourhood(i, j, sizeMatrix));
	return neighbourhood;
}

void AProceduralGenerationActor::ErosionOptAlgh( TArray<FVector> ArrayVert, float talus, TArray<FVector>& ArrayVertOut) {
	float coef = 0.5f;
	int sizeMatrix = sqrt(ArrayVert.Num());
	int i, j, k;
	int cur;
	FVector2D dmax;
	FVector2D di;
	TArray<int> neighbourhood;

	for (i = 0; i < sizeMatrix; i++)
	{
		for (j = 0; j < sizeMatrix; j++)
		{
			cur = i * sizeMatrix + j;// current cell
			dmax.X = -1;
			dmax.Y = 0;

			//neighbourhood = MoreNeighbourhood(i, j, sizeMatrix);
			//neighbourhood = ModVonNeimanNeighbourhood(i, j, sizeMatrix);
			neighbourhood = VonNeimanNeighbourhood(i, j, sizeMatrix);

			for (k = 0; k < neighbourhood.Num(); k++) {
				di.X = neighbourhood[k];
				di.Y = ArrayVert[cur].Z - ArrayVert[di.X].Z;
				if (di.Y > dmax.Y) {
					dmax.X = di.X;
					dmax.Y = di.Y;
				}
			}

			if ((dmax.X != -1)&&( 0 < dmax.Y) && (dmax.Y < talus)) {
				ArrayVert[cur].Z = ArrayVert[cur].Z -  dmax.Y/2;
				ArrayVert[dmax.X].Z += dmax.Y/2;
			}
		}
	}
	ArrayVertOut = ArrayVert;
}

TArray<FVector> AProceduralGenerationActor::voronoiDiagrams(float c1, float c2, int sizeMatrix, bool EuclideanDistanceToUse, int countRandPoint, bool domainMapUse, float parametrMinus)
{
	TArray<FVector> VertexArray;
	TArray<FVector> randomPoints;
	randomPoints.Empty();
	int i, j, k;
	FVector currentPoint;
	//create random points
	for (i = 0; i < countRandPoint - 4; i++)
	{
		currentPoint.X = FMath::RandRange(0, sizeMatrix - 1);//normal vector int to float 
		currentPoint.X = currentPoint.X / sizeMatrix;
		currentPoint.Y = FMath::RandRange(0, sizeMatrix - 1);
		currentPoint.Y = currentPoint.Y / sizeMatrix;
		currentPoint.Z = 0;
		randomPoints.Add(currentPoint);
	}
	float a = (sizeMatrix - 1);
	a = a / sizeMatrix;
	currentPoint.Set(0, 0, 0);
	randomPoints.Add(currentPoint);
	currentPoint.Set(a, 0, 0);
	randomPoints.Add(currentPoint);
	currentPoint.Set(0, a, 0);
	randomPoints.Add(currentPoint);
	currentPoint.Set(a, a, 0);
	randomPoints.Add(currentPoint);

	//generate coordinate Vert
	float minDist1;
	float minDist2;
	float dist;//distance
	int minNum1 = 0, minNum2 = 0;
	float dx, dy, d;





	for (i = 0; i < sizeMatrix; i++) {
		for (j = 0; j < sizeMatrix; j++) {
			currentPoint.X = i;
			currentPoint.X = currentPoint.X / sizeMatrix;
			currentPoint.Y = j;
			currentPoint.Y = currentPoint.Y / sizeMatrix;
			currentPoint.Z = 0;

			//find nearest points
			minDist1 = sqrt(2);
			minDist2 = minDist1;

			for (k = 0; k < countRandPoint; k++) {
				dist = FVector::Distance(randomPoints[k], currentPoint);
				if (minDist1 > dist) {
					minDist2 = minDist1;
					minNum2 = minNum1;
					minDist1 = dist;
					minNum1 = k;

				}
				else {
					if (minDist2 > dist) {
						minDist2 = dist;
						minNum2 = k;
					}
				}
			}

			//calculate h
			//calc d1
			dx = abs(currentPoint.X - randomPoints[minNum1].X);
			dy = abs(currentPoint.Y - randomPoints[minNum1].Y);
			if (EuclideanDistanceToUse)
				d = sqrt(dx * dx + dy * dy);
			else
				d = dx * dx + dy * dy;

			currentPoint.Z = c1 * d;
			//calc d2
			dx = abs(currentPoint.X - randomPoints[minNum2].X);
			dy = abs(currentPoint.Y - randomPoints[minNum2].Y);
			if (EuclideanDistanceToUse)
				d = sqrt(dx * dx + dy * dy);
			else
				d = dx * dx + dy * dy;
			currentPoint.Z += c2 * d;

			//currentPoint.Z += currentPoint.Z,;
			//currentPoint.X = ;
			//currentPoint.Y = minDist2;

			VertexArray.Add(currentPoint);

		}
	}

	//domain map
	TArray<float> domainMap;
	TArray<float> coefRandomPoint;
	for (k = 0; k < countRandPoint; k++) {
		coefRandomPoint.Add((1 - VertexArray[randomPoints[k].X * sizeMatrix + randomPoints[k].Y].Z) / 2);
	}

	for (i = 0; i < VertexArray.Num(); i++) {
		//find nearest points
		minDist1 = sqrt(2);
		for (k = 0; k < countRandPoint; k++) {
			dist = FVector::Dist2D(randomPoints[k], VertexArray[i]);
			if (minDist1 > dist) {
				minNum1 = k;
				minDist1 = dist;
			}
		}
		domainMap.Add(coefRandomPoint[minNum1]);
		VertexArray[i].Z *= coefRandomPoint[minNum1];
	}



	//road betwen halls
	for (i = 0; i < VertexArray.Num(); i++) {
		VertexArray[i].Z = VertexArray[i].Z - parametrMinus;
		if (VertexArray[i].Z < 0)
			VertexArray[i].Z = 0;
	}


	return VertexArray;
}