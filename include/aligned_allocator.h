#ifndef ALIGNED_ALLOCATOR_H
#define ALIGNED_ALLOCATOR_H

#include <cstdlib>
#include <cstdint>
#include <limits>
#include <iostream>
#include <memory>
#include <new>

/**
 * @brief Aligned allocator for SIMD vectorization.
 * 
 * Allocates memory aligned to Alignment bytes (must be power of two).
 * Uses posix_memalign on Linux/macOS, aligned_alloc on C++17, or falls back
 * to standard malloc with warning.
 */
template <typename T, size_t Alignment = 64>
class AlignedAllocator {
    static_assert(Alignment > 0 && (Alignment & (Alignment - 1)) == 0,
                  "Alignment must be a power of two");
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    
    template <typename U>
    struct rebind {
        using other = AlignedAllocator<U, Alignment>;
    };
    
    AlignedAllocator() noexcept = default;
    
    template <typename U>
    AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept {}
    
    pointer allocate(size_type n) {
        if (n > std::numeric_limits<size_type>::max() / sizeof(T)) {
            throw std::bad_alloc();
        }
        
        size_type bytes = n * sizeof(T);
        void* ptr = nullptr;
        
#if defined(_ISOC11_SOURCE) || defined(__APPLE__) || defined(__linux__)
        // Use posix_memalign on POSIX systems
        if (posix_memalign(&ptr, Alignment, bytes) != 0) {
            throw std::bad_alloc();
        }
#elif defined(_MSC_VER) || defined(__MINGW32__)
        ptr = _aligned_malloc(bytes, Alignment);
        if (!ptr) throw std::bad_alloc();
#else
        // Fallback to standard malloc (may cause SIGSEGV on SIMD loads)
        ptr = std::malloc(bytes);
        if (!ptr) throw std::bad_alloc();
        // Warn about misalignment (only in debug builds)
        #ifndef NDEBUG
        if (reinterpret_cast<uintptr_t>(ptr) % Alignment != 0) {
            std::cerr << "[WARNING] AlignedAllocator fallback may produce misaligned memory!" << std::endl;
        }
        #endif
#endif
        return static_cast<pointer>(ptr);
    }
    
    void deallocate(pointer p, size_type) noexcept {
#if defined(_ISOC11_SOURCE) || defined(__APPLE__) || defined(__linux__)
        free(p);  // posix_memalign uses free()
#elif defined(_MSC_VER) || defined(__MINGW32__)
        _aligned_free(p);
#else
        std::free(p);
#endif
    }
    
    template <typename U>
    bool operator==(const AlignedAllocator<U, Alignment>&) const noexcept {
        return true;
    }
    
    template <typename U>
    bool operator!=(const AlignedAllocator<U, Alignment>&) const noexcept {
        return false;
    }
};

#endif // ALIGNED_ALLOCATOR_H