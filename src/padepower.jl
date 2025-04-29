using LinearAlgebra
using GenericLinearAlgebra

setprecision( 256 )


function pade!( α :: AbstractVector{T}, β :: AbstractVector{T}, γ :: AbstractVector{ T }; tol = 1.0e-8 ) where T<:Real
    J = length( γ ) 
    eachindex( γ ) == 1:J || error( "The argument vector should have 1 through J indexing" )
    isodd( J ) || error( "The argument vector should have an odd number of elements" )
    K = div( J + 1, 2 )
    eachindex( β ) == eachindex( α ) == 1:K || error( "Incorrect number of elements in α or β" )
    β[1] = one( T )
    A = zeros( T, J, J )
    for c ∈ 2:K, r ∈ c:J
        A[r,c-1] = - γ[r+1-c]
    end
    for r ∈ 1:K
        A[ r, K + r - 1 ] = one( T )
    end
    # display( A )
    # λmin = GenericLinearAlgebra.eigvals( Symmetric( A'A ) )[1]
    # println( λmin )
    # λmin > tol || error( "Matrix A is not invertible (λmin = $λmin); tolerance is $tol." )
    # δ = A \ γ
    H = GenericLinearAlgebra.svd( A )
    δ = inv( A ) * γ
    β[2:K] = δ[1:K-1]
    α .= δ[K:end]
    return α, β
end

function pade( γ :: AbstractVector{ T }; tol = 1.0e-8 ) where T<:Real
    J = length( γ )
    eachindex( γ ) == 1:J || error( "The argument vector should have 1 through J indexing" )
    isodd( J ) || error( "The argument vector should have an odd number of elements" )
    K = div( J + 1, 2 )
    return pade!( zeros( T, K ), zeros( T, K), γ )
end



function ComputeDerivative( f, j, x )
    j == 0 && return f( x )
    return ForwardDiff.derivative( x->ComputeDerivative( f, j - 1, x ), x )
end


function padecoefficients( f :: Function, J :: Integer )
    return [ ComputeDerivative( f, j, 0.0 ) / factorial( j ) for j ∈ 0:J ]
end



function pade( f :: Function, J :: Integer )
    @time γ = padecoefficients( f, J )
    return pade( γ )
end




function ComputeDerivativealt( f, j, x )
    function g( x )
        ComputeDerivativealt( f, j-1, x )
    end
    return j == 0 ? f( x ) : ForwardDiff.derivative( g, x )
end


function padecoefficientsalt( f :: Function, J :: Integer )
    return [ ComputeDerivativealt( f, j, 0.0 ) / factorial( j ) for j ∈ 0:J ]
end



function padealt( f :: Function, J :: Integer )
    @time γ = padecoefficientsalt( f, J )
    println( "running pade now ")
    return pade( γ )
end


export padealt