#ifndef PARTICLEFORCES_H_
#define PARTICLEFORCES_H_
#include <iostream>
#include "typedefs.h"
#include "DynamicFactory.h"
#include "Fluid.h"

struct ParticleForceData {
	Vector fluidVelocity;
	Vector particleVelocity;
	Matrix shear;
	Fluid fluid;
	scalar particleRadius;
	scalar particleDensity;
	Vector position;
	scalar rpm;
	Vector fluidVorticity;
};

struct ParticleForce {
	virtual Vector getParticleForce(const ParticleForceData & forceData) = 0;
	virtual ParticleForce * clone() const = 0;
	virtual void writeBinary(std::ostream &) const;
	virtual void readBinary(std::istream &) { };
	virtual void fromJSON(const json &) { };
	virtual int typeId() const = 0;
};

// Particle force factory
using ParticleForceFactory = DynamicFactory<ParticleForce>;

ParticleForceFactory & particleForceFactory();

template<class DerivedType>
int registerParticleForceTypeToFactory(const std::string & typeName)
{
	return particleForceFactory().registerCreator(new NamedSimpleObjectCreator<ParticleForce, DerivedType>(typeName));
}

class StokesDrag : public ParticleForce {
public:
	Vector getParticleForce(const ParticleForceData & forceData) override;
	int typeId() const override { return typeId_; }
	StokesDrag * clone() const { return new StokesDrag(*this); }

private:
	static int typeId_;
};

class SaffmanForce : public ParticleForce {
public:
	Vector getParticleForce(const ParticleForceData & forceData) override;
	int typeId() const override { return typeId_; }
	SaffmanForce * clone() const { return new SaffmanForce(*this); }

private:
	static int typeId_;
};

class CentrifugalForce : public ParticleForce {
public:
	Vector getParticleForce(const ParticleForceData & forceData) override;
	int typeId() const override { return typeId_; }
	CentrifugalForce * clone() const { return new CentrifugalForce(*this); }

private:
	static int typeId_;
};

#endif /* PARTICLEFORCES_H_ */
