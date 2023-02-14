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
	// std::cout << (forceData.fluidVelocity - forceData.particleVelocity).norm() << std::endl;
	const scalar Re = 2. * forceData.particleRadius * (forceData.fluidVelocity - forceData.particleVelocity).norm() * forceData.particleDensity / forceData.fluid.mu();
	// std::cout << "Re: " << Re << std::endl;
	// if (Re > 1e4) {
	// 	return forceData.fluidVelocity - forceData.fluidVelocity;
	// }
	scalar Cd;
	if (Re < 0.1) {
		Cd = 300;
	} else if (Re < 1.) {
		Cd = 24. / Re;
	} else if (Re < 1000.) {
		Cd = 24. / Re * (1 + 0.15 * pow(Re, 0.687));
	} else {
		Cd = 0.44;
	}
	// std:: cout << "Cd: " << Cd << std::endl;
	Vector force = 0.5 * Cd * M_PI * forceData.fluid.rho() * forceData.particleRadius * forceData.particleRadius * (forceData.fluidVelocity - forceData.particleVelocity) * (forceData.fluidVelocity - forceData.particleVelocity).norm();
	// if (force.norm() > 1.e-11)
	// {
	// 	force.normalize();
	// 	force *= 1e-11;
	// }
	return force;
		// return 6. * M_PI * forceData.fluid.mu() * forceData.particleRadius * (forceData.fluidVelocity - forceData.particleVelocity);
	
}

int StokesDrag::typeId_ = registerParticleForceTypeToFactory<StokesDrag>("StokesDrag");

//Saffman force - stub, unusable now
Vector SaffmanForce::getParticleForce(const ParticleForceData & forceData)
{
	Vector force;
	if (forceData.fluidVorticity.norm() < 70000) {
		force = 1.61 * 4 * forceData.particleRadius * forceData.particleRadius * sqrt(forceData.fluid.mu() * forceData.fluid.rho()) * (forceData.fluidVelocity - forceData.particleVelocity).cross(forceData.fluidVorticity) / sqrt(forceData.fluidVorticity.norm());
	} else {
		const Vector vort = forceData.fluidVorticity / forceData.fluidVorticity.norm() * 30000;
		force =  1.61 * 4 * forceData.particleRadius * forceData.particleRadius * sqrt(forceData.fluid.mu() * forceData.fluid.rho()) * (forceData.fluidVelocity - forceData.particleVelocity).cross(vort) / sqrt(vort.norm());
	}
	// if (force.norm() > 1.e-11)
	// {
	// 	force.normalize();
	// 	force *= 1e-11;
	// }
	return force;
}

int SaffmanForce::typeId_ = registerParticleForceTypeToFactory<SaffmanForce>("SaffmanForce");

// Vector CentrifugalForce::getParticleForce(const ParticleForceData & forceData)
// {
// 	return 2 * M_PI * forceData.rpm/60*sqrt(forceData.position[0]*forceData.position[0] + forceData.position[1]*forceData.position[1]);
// }

// int SaffmanForce::typeId_ = registerParticleForceTypeToFactory<SaffmanForce>("SaffmanForce");
