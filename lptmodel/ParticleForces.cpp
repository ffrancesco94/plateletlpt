#include "ParticleForces.h"
#include <ostream>
#include "typedefs.h"
#include <algorithm>
#include "DynamicFactory.hh"

// Explicit instantiation of the factory
template class DynamicFactory<ParticleForce>;

ParticleForceFactory & particleForceFactory()
{
	static ParticleForceFactory particleForceFactory;
	return particleForceFactory;
}

// ParticleForce
void ParticleForce::writeBinary(std::ostream & os) const
{
	write_to_stream(os, this->typeId());
}

// StokesDrag
Vector StokesDrag::getParticleForce(const ParticleForceData & forceData)
{
	return 6. * M_PI * forceData.fluid.mu() * forceData.particleRadius * (forceData.fluidVelocity - forceData.particleVelocity);
}

int StokesDrag::typeId_ = registerParticleForceTypeToFactory<StokesDrag>("StokesDrag");

//Saffman force - stub, unusable now
Vector SaffmanForce::getParticleForce(const ParticleForceData & forceData)
{
	return 1.61 * forceData.particleRadius * forceData.particleRadius * sqrt(forceData.fluid.mu() * forceData.fluid.rho()) * (forceData.fluidVelocity - forceData.particleVelocity).cross(forceData.fluidVorticity) / sqrt(forceData.fluidVorticity.norm());
}

int SaffmanForce::typeId_ = registerParticleForceTypeToFactory<SaffmanForce>("SaffmanForce");

Vector CentrifugalForce::getParticleForce(const ParticleForceData & forceData)
{
	return 2 * M_PI * forceData.rpm/60*sqrt(forceData.position[0]*forceData.position[0] + forceData.position[1]*forceData.position[1]);
}

int SaffmanForce::typeId_ = registerParticleForceTypeToFactory<SaffmanForce>("SaffmanForce");
