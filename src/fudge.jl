


using Symbolics, PadeDiagonals, LinearAlgebra, GenericLinearAlgebra


import Base: log

setprecision( 256 )

function (@main)( args )
    @variables x 
    J = 28
    function g( x )
       log( BigFloat(1.0) + x )
    end
    f = Vector{Num}(undef, J + 1 )
    f[1] = g( x ) #atan( x )
    D = Differential( x )
    γ = zeros( BigFloat, J + 1)
    for j ∈ 2:J+1 
        f[j] = expand_derivatives( D( f[j-1 ] ), true )
        println( f[j] )
        γ[j] = build_function( f[j], x ; expression = Val{false} )( big( 0.0) ) / BigFloat( factorial( big( j - 1 ) ) )
    end
    println( γ )
    α, β = pade( γ; tol = 1.0e-32 )
    printstyled( sum( α ) / sum( β ), "   "; color = :red ); printstyled( log( big( 2.0) ), "\n"; color = :green )
    println( sum(α) / sum( β ) - log(big(2.0)) )
end
    # t = big( 0.0 )
    # g = [ build_function( f[j], x ; expression = Val{false} )( t ) for j ∈ 1:J ]
    # for j ∈ 1:J 
    #     println( g( 0.0)  / factorial( j- 1 ) )
    # end 
    
    # for j ∈ 1:J 
    #     println( build_function( f[j], x ; expression = Val{false} )( big( 0.0) ) )
    #     println( γ[j] )
    # end 
    # α, β = pade( γ; tol = 1.0e-32 )
    # println( α )
    # println( "---- " )
    # println( β )
    # printstyled( sum( α ) / sum( β ), "   "; color = :red ); printstyled( log( big( 2.0) ), "\n"; color = :green )
    # println( sum(α) / sum( β ) - log(big(2.0)) )
# end

