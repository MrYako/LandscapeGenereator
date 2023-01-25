// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "ProceduralGenerationActor.generated.h"

UCLASS()
class STRATEGYGAME_API AProceduralGenerationActor : public AActor
{
	GENERATED_BODY()
public:
			int size;
			float roughness;
			TArray<FVector> coordVert;
			TArray<float> heighmap;
public:	
	// Sets default values for this actor's properties
	AProceduralGenerationActor();
	UFUNCTION(BlueprintCallable, Category = "ProceduralTerrain", meta = (Keywords = "Voronoi Diagrams"))
		TArray<FVector> voronoiDiagrams(float c1, float c2, int sizeMatrix, bool EuclideanDistanceToUse, int countRandPoint, bool domainMapUse, float parametrMinus);
	UFUNCTION(BlueprintCallable, Category = "ProceduralTerrain", meta = (Keywords = "Diamond-Square"))
		TArray<FVector> DiamondSquareAlgh(int sizeMatrix, float roughnessCoef, bool powHeigh);
	UFUNCTION(BlueprintCallable, Category = "ProceduralTerrain", meta = (Keywords = "Combination Alghoritms"))
		TArray<FVector> CombinationAlgh(float coef1, const TArray<FVector>& coorVector1, float coef2, const TArray<FVector>& coordVector2, int sizeMatrix);
	UFUNCTION(BlueprintCallable, Category = "ProceduralTerrain", meta = (Keywords = "Thermal Erosion"))
		void thermalErosion(float talus,  TArray<FVector> ArrayVert, TArray<FVector>& ArrayVertOut);
	UFUNCTION(BlueprintCallable, Category = "ProceduralTerrain", meta = (Keywords = "Thermal Erosion calc only min"))
		void onlyMinThermalErosion(float talus, TArray<FVector> ArrayVert, TArray<FVector>& ArrayVertOut);
	UFUNCTION(BlueprintCallable, Category = "ProceduralTerrain", meta = (Keywords = "Hydraulic Erosion"))
		void hydraulicErosion( TArray<FVector> ArrayVert, TArray<FVector>& ArrayVertOut, int CountIter, TArray<FString> &Data, bool OptimizationVerson, bool WriteData);
	UFUNCTION(BlueprintCallable, Category = "ProceduralTerrain", meta = (Keywords = "Erosion Opt Algh"))
		void ErosionOptAlgh(TArray<FVector> ArrayVert, float talus, TArray<FVector>& ArrayVertOut);



	float erosionScore(const TArray<FVector>& ArrayVert);
	void DiamondSquare(int lx, int ly, int rx, int ry);
	void Diamond(int tgx, int tgy, int l);
	void Square(int lx, int ly, int rx, int ry);
	TArray<int> VonNeimanNeighbourhood(int i, int j, int sizeMatrix);
	TArray<int> ModVonNeimanNeighbourhood(int i, int j, int sizeMatrix);
	TArray<int> MoreNeighbourhood(int i, int j, int sizeMatrix);
	void transportWater(TArray<float>& material, TArray<float>& water, const TArray<FVector>& ArrayVert);
	void onlyMinTransportWater(TArray<float>& material, TArray<float>& water, const TArray<FVector>& ArrayVert);
protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

};
